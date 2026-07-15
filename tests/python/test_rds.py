from __future__ import annotations

import numpy as np

from maliampi_tools.combine import combine_seqtabs


def test_rds2py_writes_integer_matrix_with_dimnames(tmp_path):
    from rds2py import read_rds, write_rds
    from rds2py.read_matrix import MatrixWrapper

    source = tmp_path / "source.rds"
    write_rds(
        MatrixWrapper(
            np.array([[2, 0], [3, 5]], dtype=np.int32),
            dimnames=(["sample-1", "sample-2"], ["ACGT", "ACGTR"]),
        ),
        str(source),
    )
    output_h5ad = tmp_path / "pre-chimera.h5ad"
    output_rds = tmp_path / "combined.rds"
    combine_seqtabs(
        [source],
        output_h5ad,
        output_rds,
        project_id="project-a",
        dataset_id="dataset-a",
    )
    result = read_rds(str(output_rds))
    assert result.matrix.dtype == np.int32
    assert [list(names) for names in result.dimnames] == [
        ["sample-1", "sample-2"],
        ["ACGT", "ACGTR"],
    ]
