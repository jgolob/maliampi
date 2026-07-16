from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from maliampi_tools.goods import filter_goods, write_goods_outputs
from maliampi_tools.sv import SvValidationError, build_sv_artifact, read_sv_h5ad, write_sv_h5ad


def _write_source(tmp_path):
    observations = pd.DataFrame(
        {
            "project_id": ["project"] * 3,
            "deposited_specimen_id": ["A", "B", "C"],
        }
    )
    artifact = build_sv_artifact(
        sparse.csr_matrix([[12, 4, 0], [10, 4, 0], [1, 0, 1]], dtype=np.int64),
        observations,
        ["ACGT", "ACGA", "ACGC"],
        dataset_id="goods-test",
    )
    source = tmp_path / "sv.h5ad"
    write_sv_h5ad(artifact, source)
    return source


def test_goods_is_deterministic_and_records_outputs(tmp_path):
    source = _write_source(tmp_path)
    parameters = {
        "convergence_delta": 1.0,
        "min_reads": 1,
        "min_prevalence": 1,
        "seed": 17,
    }
    first, first_convergence, first_curves = filter_goods(source, **parameters)
    second, second_convergence, second_curves = filter_goods(source, **parameters)

    np.testing.assert_array_equal(first.X.toarray(), second.X.toarray())
    pd.testing.assert_frame_equal(first_convergence, second_convergence)
    pd.testing.assert_frame_equal(first_curves, second_curves)
    assert first.uns["maliampi"]["provenance"]["goods_filter"]["seed"] == 17
    assert set(first.obs["goods_disposition"]) == {"converged"}
    assert first_curves.groupby("observation_id").size().sum() == 32

    output = tmp_path / "sv.goods.h5ad"
    convergence = tmp_path / "goods.convergence.parquet"
    curves = tmp_path / "goods.curves.parquet"
    write_goods_outputs(source, output, convergence, curves, **parameters)
    read_sv_h5ad(output)
    assert len(pd.read_parquet(convergence)) == 3
    assert not pd.read_parquet(curves).empty


def test_goods_nonconverged_specimens_do_not_contribute_to_prevalence(tmp_path):
    source = _write_source(tmp_path)
    dropped, disposition, _ = filter_goods(
        source,
        convergence_delta=1.0,
        min_reads=3,
        min_prevalence=1,
        keep_nonconverged=False,
    )
    kept, _, _ = filter_goods(
        source,
        convergence_delta=1.0,
        min_reads=3,
        min_prevalence=1,
        keep_nonconverged=True,
    )

    assert dropped.n_obs == 2
    assert kept.n_obs == 3
    assert not disposition.set_index("deposited_specimen_id").loc["C", "goods_converged"]
    assert kept.obs.set_index("deposited_specimen_id").loc["C", "goods_disposition"] == (
        "retained_nonconverged"
    )
    assert "ACGC" not in kept.var["sequence"].tolist()


@pytest.mark.parametrize(
    ("parameter", "value"),
    [("convergence_delta", -0.1), ("min_reads", 0), ("min_prevalence", 0)],
)
def test_goods_rejects_invalid_cutoffs(tmp_path, parameter, value):
    source = _write_source(tmp_path)
    with pytest.raises(SvValidationError):
        filter_goods(source, **{parameter: value})


def test_goods_rejects_empty_result(tmp_path):
    source = _write_source(tmp_path)
    with pytest.raises(SvValidationError, match="removed every"):
        filter_goods(
            source,
            convergence_delta=1.0,
            min_reads=1,
            min_prevalence=4,
        )
