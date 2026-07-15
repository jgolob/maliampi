"""Creation and validation of the versioned sparse SV H5AD contract."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from hashlib import sha256
from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SV_SCHEMA_VERSION = "1.0.0"
IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


class SvValidationError(ValueError):
    """Raised when an SV artifact does not meet the MaliAmPi contract."""


def normalize_sequence(sequence: str) -> str:
    """Normalize a sequence and reject characters outside uppercase IUPAC DNA."""
    normalized = "".join(str(sequence).split()).upper()
    if not normalized or any(base not in IUPAC_DNA for base in normalized):
        raise SvValidationError(f"Invalid IUPAC DNA sequence: {sequence!r}")
    return normalized


def sequence_hash(sequence: str) -> str:
    return sha256(normalize_sequence(sequence).encode("ascii")).hexdigest()


def sv_id_from_hash(digest: str) -> str:
    """Create a readable, deterministic ID from disjoint SHA-256 portions."""
    if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
        raise SvValidationError("SV SHA-256 digest must be 64 lowercase hexadecimal characters")
    decimal = int(digest[:8], 16)
    suffix = digest[8:24]
    return f"sv-{decimal:010d}-{suffix}"


def observation_id(project_id: str, specimen_id: str) -> tuple[str, str]:
    """Return the persistent observation name and full identity digest."""
    project, specimen = str(project_id), str(specimen_id)
    if not project or not specimen:
        raise SvValidationError("project_id and deposited specimen_id must not be empty")
    normalized = f"{project.strip()}\0{specimen.strip()}".encode()
    digest = sha256(normalized).hexdigest()
    return f"{project}/{specimen}~{digest[:16]}", digest


def _as_csr_counts(matrix: Any) -> sparse.csr_matrix:
    csr = sparse.csr_matrix(matrix)
    if not np.isfinite(csr.data).all() or np.any(csr.data < 0):
        raise SvValidationError("SV counts must be finite and nonnegative")
    if not np.equal(csr.data, np.floor(csr.data)).all():
        raise SvValidationError("SV counts must be integral")
    if csr.data.size and csr.data.max() > np.iinfo(np.int64).max:
        raise SvValidationError("SV count exceeds int64 range")
    csr = csr.astype(np.int64)
    csr.eliminate_zeros()
    return csr


def build_sv_artifact(
    counts: Any,
    observations: pd.DataFrame,
    sequences: Iterable[str],
    *,
    dataset_id: str,
    tool_version: str = "unknown",
    provenance: Mapping[str, Any] | None = None,
    refpkg_identity: str | None = None,
) -> ad.AnnData:
    """Build a validated specimen-by-SV artifact from sparse-compatible inputs."""
    sequence_list = [normalize_sequence(sequence) for sequence in sequences]
    hashes = [sequence_hash(sequence) for sequence in sequence_list]
    sv_ids = [sv_id_from_hash(digest) for digest in hashes]
    if len(set(hashes)) != len(hashes) or len(set(sv_ids)) != len(sv_ids):
        raise SvValidationError("SV sequences and composite IDs must be unique")
    required = {"project_id", "deposited_specimen_id"}
    if not required.issubset(observations.columns):
        raise SvValidationError(f"observations must contain {sorted(required)}")
    obs = observations.copy().reset_index(drop=True)
    identities = [
        observation_id(project, specimen)
        for project, specimen in zip(obs["project_id"], obs["deposited_specimen_id"], strict=True)
    ]
    obs.index = pd.Index([identity[0] for identity in identities], name="observation_id")
    obs["identity_sha256"] = [identity[1] for identity in identities]
    if not obs.index.is_unique:
        raise SvValidationError("A project/specimen identity occurs more than once")
    matrix = _as_csr_counts(counts)
    if matrix.shape != (len(obs), len(sequence_list)):
        raise SvValidationError("Count matrix dimensions do not match observations and sequences")
    var = pd.DataFrame(
        {
            "sv_id": sv_ids,
            "sequence": sequence_list,
            "sequence_sha256": hashes,
            "sequence_length": [len(sequence) for sequence in sequence_list],
            "total_abundance": np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64),
        },
        index=pd.Index(sv_ids, name="sv_id"),
    )
    metadata: dict[str, Any] = {
        "schema_version": SV_SCHEMA_VERSION,
        "dataset_id": str(dataset_id),
        "tool_version": str(tool_version),
        "provenance": dict(provenance or {}),
    }
    if refpkg_identity is not None:
        metadata["refpkg_identity"] = str(refpkg_identity)
    artifact = ad.AnnData(X=matrix, obs=obs, var=var, uns={"maliampi": metadata})
    validate_sv_artifact(artifact)
    return artifact


def validate_sv_artifact(artifact: ad.AnnData) -> None:
    """Validate all invariants required whenever an SV H5AD is read."""
    if not sparse.isspmatrix_csr(artifact.X):
        raise SvValidationError("sv.h5ad X must be CSR")
    matrix = _as_csr_counts(artifact.X)
    if artifact.n_obs == 0 or artifact.n_vars == 0:
        raise SvValidationError("sv.h5ad must not contain orphan specimens or SVs")
    if np.any(np.asarray(matrix.sum(axis=1)).ravel() == 0):
        raise SvValidationError("sv.h5ad contains an orphan specimen with no SV counts")
    if np.any(np.asarray(matrix.sum(axis=0)).ravel() == 0):
        raise SvValidationError("sv.h5ad contains an orphan SV with no counts")
    required_obs = {"project_id", "deposited_specimen_id", "identity_sha256"}
    required_var = {"sv_id", "sequence", "sequence_sha256", "sequence_length", "total_abundance"}
    if not required_obs.issubset(artifact.obs.columns):
        raise SvValidationError("sv.h5ad obs lacks required specimen identity fields")
    if not required_var.issubset(artifact.var.columns):
        raise SvValidationError("sv.h5ad var lacks required SV fields")
    if not artifact.obs_names.is_unique or not artifact.var_names.is_unique:
        raise SvValidationError("sv.h5ad observation and variable names must be unique")
    obs = cast(pd.DataFrame, artifact.obs)
    var = cast(pd.DataFrame, artifact.var)
    expected_obs = [
        observation_id(project, specimen)
        for project, specimen in zip(obs["project_id"], obs["deposited_specimen_id"], strict=True)
    ]
    if artifact.obs_names.tolist() != [item[0] for item in expected_obs]:
        raise SvValidationError("sv.h5ad observation names do not match deposited identity")
    if obs["identity_sha256"].tolist() != [item[1] for item in expected_obs]:
        raise SvValidationError("sv.h5ad observation hashes do not match deposited identity")
    computed_hashes = [sequence_hash(sequence) for sequence in var["sequence"]]
    if var["sequence_sha256"].tolist() != computed_hashes:
        raise SvValidationError("sv.h5ad sequence hashes do not match sequences")
    expected_ids = [sv_id_from_hash(digest) for digest in computed_hashes]
    if var["sv_id"].tolist() != expected_ids or artifact.var_names.tolist() != expected_ids:
        raise SvValidationError("sv.h5ad composite SV IDs do not match sequence hashes")
    if var["sequence_length"].tolist() != [len(sequence) for sequence in var["sequence"]]:
        raise SvValidationError("sv.h5ad sequence lengths do not match sequences")
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64).tolist()
    if var["total_abundance"].astype(np.int64).tolist() != totals:
        raise SvValidationError("sv.h5ad total abundances do not match X")
    metadata = artifact.uns.get("maliampi")
    if not isinstance(metadata, dict) or metadata.get("schema_version") != SV_SCHEMA_VERSION:
        raise SvValidationError("sv.h5ad lacks a supported maliampi schema declaration")
    if not metadata.get("dataset_id") or not metadata.get("tool_version"):
        raise SvValidationError("sv.h5ad lacks dataset or tool identity")


def read_sv_h5ad(path: str | Path) -> ad.AnnData:
    artifact = ad.read_h5ad(path)
    validate_sv_artifact(artifact)
    return artifact


def write_sv_h5ad(artifact: ad.AnnData, path: str | Path) -> None:
    validate_sv_artifact(artifact)
    cast(Any, artifact).write_h5ad(str(path))


def filter_sv_artifact(artifact: ad.AnnData, retained_hashes: Iterable[str]) -> ad.AnnData:
    """Filter a pre-chimera artifact by retained full sequence hashes."""
    retained = set(retained_hashes)
    missing = retained - set(artifact.var["sequence_sha256"])
    if missing:
        raise SvValidationError("DADA2 retained a sequence not present in pre-chimera sv.h5ad")
    filtered = artifact[:, artifact.var["sequence_sha256"].isin(retained)].copy()
    matrix = _as_csr_counts(filtered.X)
    filtered.var["total_abundance"] = np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64)
    filtered.uns["maliampi"]["provenance"] = dict(filtered.uns["maliampi"].get("provenance", {}))
    filtered.uns["maliampi"]["provenance"]["chimera_filter"] = "applied"
    validate_sv_artifact(filtered)
    return filtered
