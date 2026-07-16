from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd
import pytest

from maliampi_tools.sv import SvValidationError, read_sv_h5ad
from maliampi_tools.swarm import (
    combine_dsv_registries,
    filter_swarm_seeds,
    finalize_swarm_h5ad,
    prepare_specimen_dsv,
)


def _write_fasta(path, records):
    with gzip.open(path, "wt") as handle:
        for name, sequence, count in records:
            handle.write(f">{name};size={count}\n{sequence}\n")


def _prepare_registry(tmp_path, specimen, batch, records, suffix):
    derep = tmp_path / f"{suffix}.fasta.gz"
    registry = tmp_path / f"{suffix}.parquet"
    _write_fasta(derep, records)
    prepare_specimen_dsv(derep, registry, specimen=specimen, batch=batch)
    return registry


def test_swarm_builds_sparse_h5ad_and_aggregates_repeated_runs(tmp_path):
    registries = [
        _prepare_registry(tmp_path, "A", "run-1", [("a", "ACGT", 3), ("b", "ACGA", 2)], "a1"),
        _prepare_registry(tmp_path, "A", "run-2", [("a", "ACGT", 4)], "a2"),
        _prepare_registry(tmp_path, "B", "run-1", [("a", "ACGC", 5)], "b1"),
    ]
    dsv_fasta = tmp_path / "dsv.fasta.gz"
    registry = tmp_path / "registry.parquet"
    combine_dsv_registries(registries, dsv_fasta, registry)
    table = pd.read_parquet(registry)
    dsv = table.drop_duplicates("sequence").set_index("sequence")["dsv_id"].to_dict()

    clusters = tmp_path / "clusters.txt"
    clusters.write_text(f"{dsv['ACGT']};size=7 {dsv['ACGA']};size=2\n{dsv['ACGC']};size=5\n")
    nonchimera = tmp_path / "nonchimera.fasta.gz"
    _write_fasta(nonchimera, [(dsv["ACGT"], "ACGT", 9), (dsv["ACGC"], "ACGC", 5)])
    output = tmp_path / "sv.h5ad"
    stats = tmp_path / "swarm.stats.parquet"
    finalize_swarm_h5ad(
        registry,
        clusters,
        nonchimera,
        output,
        stats,
        project_id="project",
        dataset_id="dataset",
    )

    artifact = read_sv_h5ad(output)
    assert artifact.obs["deposited_specimen_id"].tolist() == ["A", "B"]
    assert json.loads(artifact.obs.iloc[0]["batches"]) == ["run-1", "run-2"]
    np.testing.assert_array_equal(artifact.X.toarray(), [[9, 0], [0, 5]])
    summary = pd.read_parquet(stats).iloc[0]
    assert summary["input_reads"] == 14
    assert summary["retained_reads"] == 14


def test_swarm_rejects_malformed_cluster_coverage(tmp_path):
    registry = _prepare_registry(tmp_path, "A", "run", [("a", "ACGT", 3)], "a")
    dsv_fasta = tmp_path / "dsv.fasta.gz"
    combined = tmp_path / "combined.parquet"
    combine_dsv_registries([registry], dsv_fasta, combined)
    clusters = tmp_path / "clusters.txt"
    clusters.write_text("unknown;size=3\n")
    nonchimera = tmp_path / "nonchimera.fasta.gz"
    _write_fasta(nonchimera, [("unknown", "ACGT", 3)])

    with pytest.raises(SvValidationError, match="do not cover"):
        finalize_swarm_h5ad(
            combined,
            clusters,
            nonchimera,
            tmp_path / "sv.h5ad",
            tmp_path / "stats.parquet",
            project_id="project",
            dataset_id="dataset",
        )


def test_swarm_singleton_filter_and_empty_rejection(tmp_path):
    seeds = tmp_path / "seeds.fasta.gz"
    _write_fasta(seeds, [("keep", "ACGT", 2), ("drop", "ACGA", 1)])
    output = tmp_path / "filtered.fasta.gz"
    filter_swarm_seeds(seeds, output)
    with gzip.open(output, "rt") as handle:
        filtered = handle.read()
    assert "keep;size=2" in filtered
    assert "drop" not in filtered

    singleton = tmp_path / "singleton.fasta.gz"
    _write_fasta(singleton, [("drop", "ACGA", 1)])
    with pytest.raises(SvValidationError, match="removed every"):
        filter_swarm_seeds(singleton, tmp_path / "none.fasta.gz")
