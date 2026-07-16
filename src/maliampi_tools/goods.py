"""Deterministic Good's-coverage filtering for sparse SV artifacts."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .sv import SvValidationError, read_sv_h5ad, validate_sv_artifact, write_sv_h5ad

GOODS_ALGORITHM_VERSION = "1.0.0"


def _observation_seed(seed: int, observation_id: str) -> int:
    digest = sha256(f"{seed}\0{observation_id}".encode()).digest()
    return int.from_bytes(digest[:8], "big")


def filter_goods(
    source_path: str | Path,
    *,
    convergence_delta: float = 0.0001,
    min_reads: int = 30,
    min_prevalence: int = 2,
    seed: int = 1,
    keep_nonconverged: bool = False,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    """Filter an SV artifact using parameterized Good's-coverage convergence."""
    if not 0 <= convergence_delta <= 1:
        raise SvValidationError("Good's convergence delta must be between zero and one")
    if min_reads < 1 or min_prevalence < 1:
        raise SvValidationError("Good's minimum reads and prevalence must be positive")

    source_path = Path(source_path)
    source = read_sv_h5ad(source_path)
    matrix = sparse.csr_matrix(source.X).astype(np.int64)
    filtered = matrix.copy().tolil()
    convergence_rows: list[dict[str, object]] = []
    curve_rows: list[dict[str, object]] = []
    converged = np.zeros(source.n_obs, dtype=bool)
    threshold_reads = np.full(source.n_obs, -1, dtype=np.int64)

    for obs_index, observation_id in enumerate(source.obs_names):
        counts = matrix.getrow(obs_index).toarray().ravel().astype(np.int64)
        read_order = np.repeat(np.arange(source.n_vars, dtype=np.int64), counts)
        rng = np.random.default_rng(_observation_seed(seed, observation_id))
        rng.shuffle(read_order)
        seen = np.zeros(source.n_vars, dtype=np.int64)
        singleton_count = 0
        observed_count = 0
        previous_coverage: float | None = None
        threshold_read: int | None = None

        for read_number, sv_index in enumerate(read_order, start=1):
            old_count = seen[sv_index]
            seen[sv_index] += 1
            if old_count == 0:
                observed_count += 1
                singleton_count += 1
            elif old_count == 1:
                singleton_count -= 1
            coverage = 1.0 - singleton_count / read_number
            curve_rows.append(
                {
                    "observation_id": observation_id,
                    "read_number": read_number,
                    "sv_observed": observed_count,
                    "singleton_count": singleton_count,
                    "goods_coverage": coverage,
                }
            )
            if (
                threshold_read is None
                and previous_coverage is not None
                and read_number > min_reads
                and abs(coverage - previous_coverage) <= convergence_delta
            ):
                threshold_read = read_number
            previous_coverage = coverage

        if threshold_read is not None:
            converged[obs_index] = True
            threshold_reads[obs_index] = threshold_read
            abundance_cutoff = threshold_read / 4.0
            filtered[obs_index, :] = np.where(counts >= abundance_cutoff, counts, 0)
        convergence_rows.append(
            {
                "observation_id": observation_id,
                "project_id": source.obs.iloc[obs_index]["project_id"],
                "deposited_specimen_id": source.obs.iloc[obs_index]["deposited_specimen_id"],
                "goods_converged": bool(converged[obs_index]),
                "threshold_read": threshold_read,
                "total_reads": int(counts.sum()),
                "retained_before_prevalence": int(filtered.getrow(obs_index).count_nonzero()),
            }
        )

    source.obs["goods_converged"] = converged
    source.obs["goods_threshold_read"] = threshold_reads
    source.obs["goods_disposition"] = np.where(
        converged, "converged", "retained_nonconverged" if keep_nonconverged else "dropped"
    )
    filtered_csr = filtered.tocsr()
    prevalence = np.asarray((filtered_csr[converged, :] > 0).sum(axis=0)).ravel()
    keep_vars = prevalence >= min_prevalence
    keep_obs = converged.copy() if not keep_nonconverged else np.ones(source.n_obs, dtype=bool)
    result = source[keep_obs, keep_vars].copy()
    result_matrix = filtered_csr[keep_obs, :][:, keep_vars]
    result.X = result_matrix
    row_nonzero = np.asarray(result_matrix.sum(axis=1)).ravel() > 0
    column_nonzero = np.asarray(result_matrix.sum(axis=0)).ravel() > 0
    result = result[row_nonzero, column_nonzero].copy()
    if result.n_obs == 0 or result.n_vars == 0:
        raise SvValidationError("Good's filtering removed every specimen or SV")
    result_matrix = sparse.csr_matrix(result.X)
    result.var["total_abundance"] = np.asarray(result_matrix.sum(axis=0)).ravel().astype(np.int64)
    provenance = dict(result.uns["maliampi"].get("provenance", {}))
    provenance["goods_filter"] = {
        "algorithm_version": GOODS_ALGORITHM_VERSION,
        "source_sv_artifact_digest": sha256(source_path.read_bytes()).hexdigest(),
        "convergence_delta": convergence_delta,
        "min_reads": min_reads,
        "min_prevalence": min_prevalence,
        "seed": seed,
        "keep_nonconverged": keep_nonconverged,
        "converged_observations": int(converged.sum()),
        "nonconverged_observations": int((~converged).sum()),
        "retained_observations": int(result.n_obs),
        "dropped_observations": int(source.n_obs - result.n_obs),
        "input_observations": int(source.n_obs),
        "nonconverged_observation_ids": source.obs_names[~converged].tolist(),
        "retained_observation_ids": result.obs_names.tolist(),
    }
    result.uns["maliampi"]["provenance"] = provenance
    validate_sv_artifact(result)
    return result, pd.DataFrame(convergence_rows), pd.DataFrame(curve_rows)


def write_goods_outputs(
    source_path: str | Path,
    output_h5ad: str | Path,
    convergence_parquet: str | Path,
    curves_parquet: str | Path,
    **parameters: Any,
) -> None:
    """Filter and write all canonical Good's outputs."""
    result, convergence, curves = filter_goods(source_path, **parameters)
    write_sv_h5ad(result, output_h5ad)
    convergence.to_parquet(convergence_parquet, index=False)
    curves.to_parquet(curves_parquet, index=False)
