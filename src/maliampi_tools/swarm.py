"""Sparse registry and H5AD adapters for the external Swarm engine."""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path
from typing import TextIO, cast

import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy import sparse

from .sv import SvValidationError, build_sv_artifact, normalize_sequence, write_sv_h5ad

SIZE_PATTERN = re.compile(r";size=(\d+)(?:;|$)")


def _open_text(path: str | Path, mode: str = "rt") -> TextIO:
    handle = gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)  # noqa: SIM115
    return cast(TextIO, handle)


def prepare_specimen_dsv(
    derep_fasta: str | Path,
    output_registry: str | Path,
    *,
    specimen: str,
    batch: str,
) -> None:
    """Normalize one specimen's vsearch dereplication output."""
    rows: list[dict[str, object]] = []
    with _open_text(derep_fasta) as handle:
        for index, record in enumerate(SeqIO.parse(handle, "fasta"), start=1):
            match = SIZE_PATTERN.search(record.id)
            if match is None:
                raise SvValidationError(f"Dereplicated record lacks size annotation: {record.id}")
            sequence = normalize_sequence(str(record.seq))
            local_id = f"local-{index:010d}"
            count = int(match.group(1))
            if count <= 0:
                raise SvValidationError(f"Dereplicated record has a nonpositive size: {record.id}")
            rows.append(
                {
                    "specimen": specimen,
                    "batch": batch,
                    "local_id": local_id,
                    "sequence": sequence,
                    "count": count,
                }
            )
    if not rows:
        raise SvValidationError(f"Specimen {specimen!r} has no dereplicated sequences")
    pd.DataFrame(rows).to_parquet(output_registry, index=False)


def combine_dsv_registries(
    registries: list[str | Path], output_fasta: str | Path, output_registry: str | Path
) -> None:
    """Create deterministic global DSV identities and a sparse-compatible registry."""
    if not registries:
        raise SvValidationError("Swarm requires at least one specimen registry")
    combined = pd.concat((pd.read_parquet(path) for path in registries), ignore_index=True)
    required = {"specimen", "batch", "sequence", "count"}
    if required - set(combined.columns):
        raise SvValidationError("Swarm specimen registry lacks required columns")
    if combined[list(required)].isna().any().any():
        raise SvValidationError("Swarm specimen registry contains missing values")
    if (combined["count"] <= 0).any() or (combined["count"] % 1 != 0).any():
        raise SvValidationError("Swarm registry counts must be positive integers")
    totals = combined.groupby("sequence", as_index=False)["count"].sum()
    totals = totals.sort_values(["count", "sequence"], ascending=[False, True]).reset_index(
        drop=True
    )
    totals["dsv_id"] = [f"DSV{index:010d}" for index in range(1, len(totals) + 1)]
    dsv_for_sequence = totals.set_index("sequence")["dsv_id"]
    combined["dsv_id"] = combined["sequence"].map(dsv_for_sequence)
    combined.to_parquet(output_registry, index=False)
    with _open_text(output_fasta, "wt") as output:
        for sequence, count, dsv_id in totals[["sequence", "count", "dsv_id"]].itertuples(
            index=False, name=None
        ):
            output.write(f">{dsv_id};size={int(count)}\n{sequence}\n")


def filter_swarm_seeds(input_fasta: str | Path, output_fasta: str | Path) -> None:
    """Remove global singleton seeds before chimera detection."""
    retained = 0
    with _open_text(input_fasta) as handle, _open_text(output_fasta, "wt") as output:
        for record in SeqIO.parse(handle, "fasta"):
            match = SIZE_PATTERN.search(record.id)
            if match is None:
                raise SvValidationError(f"Swarm seed lacks size annotation: {record.id}")
            if int(match.group(1)) > 1:
                output.write(f">{record.id}\n{normalize_sequence(str(record.seq))}\n")
                retained += 1
    if retained == 0:
        raise SvValidationError("Swarm singleton filtering removed every seed")


def finalize_swarm_h5ad(
    registry_parquet: str | Path,
    clusters_path: str | Path,
    nonchimera_fasta: str | Path,
    output_h5ad: str | Path,
    stats_parquet: str | Path,
    *,
    project_id: str,
    dataset_id: str,
) -> None:
    """Aggregate DSV counts by retained Swarm clusters into canonical H5AD."""
    registry = pd.read_parquet(registry_parquet)
    required = {"specimen", "batch", "sequence", "count", "dsv_id"}
    if required - set(registry.columns):
        raise SvValidationError("Swarm combined registry lacks required columns")
    dsv_to_cluster: dict[str, int] = {}
    with open(clusters_path, encoding="utf-8") as handle:
        for cluster_index, line in enumerate(handle):
            members = [token.split(";", 1)[0] for token in line.split()]
            if not members:
                raise SvValidationError("Swarm cluster file contains an empty cluster")
            for dsv_id in members:
                if dsv_id in dsv_to_cluster:
                    raise SvValidationError(f"DSV occurs in multiple Swarm clusters: {dsv_id}")
                dsv_to_cluster[dsv_id] = cluster_index
    if set(registry["dsv_id"]) - set(dsv_to_cluster):
        raise SvValidationError("Swarm clusters do not cover every registered DSV")

    retained_sequences: dict[int, str] = {}
    with _open_text(nonchimera_fasta) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seed_id = record.id.split(";", 1)[0]
            if seed_id not in dsv_to_cluster:
                raise SvValidationError(f"Nonchimera seed is absent from clusters: {seed_id}")
            cluster = dsv_to_cluster[seed_id]
            if cluster in retained_sequences:
                raise SvValidationError("Nonchimera FASTA contains multiple seeds for one cluster")
            retained_sequences[cluster] = normalize_sequence(str(record.seq))
    if not retained_sequences:
        raise SvValidationError("Swarm chimera filtering retained no clusters")

    registry["cluster"] = registry["dsv_id"].map(dsv_to_cluster)
    kept = registry[registry["cluster"].isin(retained_sequences)].copy()
    specimen_order = sorted(pd.unique(registry["specimen"]))
    cluster_order = sorted(retained_sequences)
    specimen_index = {specimen: index for index, specimen in enumerate(specimen_order)}
    cluster_index = {cluster: index for index, cluster in enumerate(cluster_order)}
    aggregated = kept.groupby(["specimen", "cluster"], as_index=False)["count"].sum()
    counts = sparse.coo_matrix(
        (
            aggregated["count"].astype(np.int64),
            (
                aggregated["specimen"].map(specimen_index),
                aggregated["cluster"].map(cluster_index),
            ),
        ),
        shape=(len(specimen_order), len(cluster_order)),
        dtype=np.int64,
    ).tocsr()
    row_nonzero = np.asarray(counts.sum(axis=1)).ravel() > 0
    observations = (
        pd.DataFrame(
            {
                "project_id": project_id,
                "deposited_specimen_id": specimen_order,
                "batches": [
                    json.dumps(
                        sorted(
                            pd.unique(
                                registry.loc[registry["specimen"] == specimen, "batch"]
                            ).tolist()
                        )
                    )
                    for specimen in specimen_order
                ],
            }
        )
        .loc[row_nonzero]
        .reset_index(drop=True)
    )
    counts = counts[row_nonzero, :]
    artifact = build_sv_artifact(
        counts,
        observations,
        [retained_sequences[cluster] for cluster in cluster_order],
        dataset_id=dataset_id,
        tool_version="swarm",
        provenance={"sv_engine": "swarm", "cluster_count": len(cluster_order)},
    )
    write_sv_h5ad(artifact, output_h5ad)
    input_reads = int(registry["count"].sum())
    kept_reads = int(kept["count"].sum())
    pd.DataFrame(
        [
            {
                "input_dsvs": int(registry["dsv_id"].nunique()),
                "input_reads": input_reads,
                "retained_clusters": len(cluster_order),
                "retained_reads": kept_reads,
                "filtered_reads": input_reads - kept_reads,
            }
        ]
    ).to_parquet(stats_parquet, index=False)
