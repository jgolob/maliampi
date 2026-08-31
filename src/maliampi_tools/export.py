"""Transient legacy exports derived from the canonical SV artifact."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy import sparse

from .sv import read_sv_h5ad


def export_legacy(
    h5ad: str | Path, output_dir: str | Path, *, pplacer_reduplication: bool = False
) -> None:
    """Materialize compatibility files only at an external-tool boundary."""
    artifact = read_sv_h5ad(h5ad)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    sv_ids = artifact.var["sv_id"].to_numpy()
    sequences = artifact.var["sequence"].to_numpy()
    obs_names = artifact.obs_names.to_numpy()

    # FASTA — streaming, no memory overhead
    with (output / "sv.fasta").open("w", encoding="ascii") as handle:
        for sv_id, sequence in zip(sv_ids, sequences, strict=True):
            handle.write(f">{sv_id}\n{sequence}\n")

    # Build long table directly from sparse COO coordinates — no list-of-dicts
    matrix = sparse.csr_matrix(artifact.X).tocoo()
    specimens = obs_names[matrix.row]
    svs = sv_ids[matrix.col]  # type: ignore[index]
    counts = matrix.data.astype(np.int64)

    # sv.long.csv — stream to disk
    with (output / "sv.long.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("specimen", "sv_id", "abundance"))
        for specimen, sv_id, count in zip(specimens, svs, counts, strict=True):
            writer.writerow((specimen, sv_id, int(count)))

    # sv.multiplicity.csv — headerless, reordered columns for gappa
    with (output / "sv.multiplicity.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        for sv_id, specimen, count in zip(svs, specimens, counts, strict=True):
            writer.writerow((sv_id, specimen, int(count)))

    # sv.share.csv — pivot via sparse-to-dense on the original matrix (specimens x SVs)
    # Only materialize the dense matrix, not a redundant DataFrame
    with (output / "sv.share.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sv_id", *obs_names])
        dense = matrix.tocsr().toarray()  # shape: (n_obs, n_var)
        for var_idx, sv_id in enumerate(sv_ids):
            writer.writerow([sv_id, *(int(v) for v in dense[:, var_idx])])

    if pplacer_reduplication:
        # sv.weights.csv — per-SV total abundance
        totals = np.asarray(sparse.csr_matrix(artifact.X).sum(axis=0)).ravel().astype(np.int64)
        with (output / "sv.weights.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(("sv_id", "weight"))
            for sv_id, weight in zip(sv_ids, totals, strict=True):
                writer.writerow((sv_id, int(weight)))

        # sv.map.csv — headerless sv_id:specimen mapping
        with (output / "sv.map.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            for sv_id, specimen in zip(svs, specimens, strict=True):
                writer.writerow((f"{sv_id}:{specimen}",))
