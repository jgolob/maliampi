"""Transient legacy exports derived from the canonical SV artifact."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy import sparse

from .sv import read_sv_h5ad


def export_legacy(
    h5ad: str | Path, output_dir: str | Path, *, pplacer_reduplication: bool = False
) -> None:
    """Materialize compatibility files only at an external-tool boundary."""
    artifact = read_sv_h5ad(h5ad)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    var = artifact.var.copy()
    with (output / "sv.fasta").open("w", encoding="ascii") as handle:
        for sv_id, sequence in zip(var["sv_id"], var["sequence"], strict=True):
            handle.write(f">{sv_id}\n{sequence}\n")
    long_rows: list[dict[str, object]] = []
    matrix = sparse.csr_matrix(artifact.X).tocoo()
    for observation_index, sv_index, value in zip(matrix.row, matrix.col, matrix.data, strict=True):
        count = int(value)
        long_rows.append(
            {
                "specimen": artifact.obs_names[observation_index],
                "sv_id": artifact.var_names[sv_index],
                "abundance": count,
            }
        )
    long = pd.DataFrame(long_rows).reindex(columns=("specimen", "sv_id", "abundance"))
    long.to_csv(output / "sv.long.csv", index=False)
    share = (
        long.pivot(index="sv_id", columns="specimen", values="abundance").fillna(0).astype("int64")
    )
    share.to_csv(output / "sv.share.csv")
    long.loc[:, ["sv_id", "specimen", "abundance"]].to_csv(
        output / "sv.multiplicity.csv", index=False, header=False
    )
    if pplacer_reduplication:
        multiplicities = long.groupby("sv_id", as_index=False)["abundance"].sum()
        weights = multiplicities.rename(columns={"abundance": "weight"})
        weights.to_csv(output / "sv.weights.csv", index=False)
        mapping = long.assign(map=lambda frame: frame["sv_id"] + ":" + frame["specimen"])["map"]
        mapping.to_csv(output / "sv.map.csv", index=False, header=False)
