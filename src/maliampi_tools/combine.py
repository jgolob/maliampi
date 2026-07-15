"""Incremental sparse DADA2 sequence-table combination and RDS exchange."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from .sv import SvValidationError, build_sv_artifact, write_sv_h5ad


def combine_seqtabs(
    rds_paths: list[str | Path],
    output_h5ad: str | Path,
    output_rds: str | Path,
    *,
    project_id: str,
    dataset_id: str,
) -> None:
    """Read RDS tables one at a time and create CSR plus a DADA2 transient RDS."""
    try:
        from rds2py import read_rds, write_rds
        from rds2py.read_matrix import MatrixWrapper
    except ImportError as error:  # pragma: no cover - dependency is mandatory
        raise RuntimeError("rds2py is required for the DADA2 RDS compatibility bridge") from error
    sequence_index: dict[str, int] = {}
    specimen_rows: dict[str, dict[int, int]] = defaultdict(dict)
    specimen_runs: dict[str, set[str]] = {}
    for path in map(Path, rds_paths):
        value = read_rds(str(path))
        matrix = getattr(value, "matrix", value)
        dimnames = getattr(value, "dimnames", None)
        if dimnames is None or len(dimnames) != 2:
            raise SvValidationError(f"{path} is not a DADA2 sequence table with dimnames")
        samples, sequences = dimnames
        table = sparse.csr_matrix(matrix)
        if table.shape != (len(samples), len(sequences)):
            raise SvValidationError(f"{path} dimensions do not agree with RDS dimnames")
        for row, specimen in enumerate(samples):
            specimen = str(specimen)
            specimen_runs.setdefault(specimen, set()).add(str(path))
            for column, count in zip(table[row].indices, table[row].data, strict=True):
                if not float(count).is_integer() or count < 0:
                    raise SvValidationError(f"{path} has non-integral or negative counts")
                sequence = str(sequences[column])
                global_column = sequence_index.setdefault(sequence, len(sequence_index))
                specimen_rows[specimen][global_column] = specimen_rows[specimen].get(
                    global_column, 0
                ) + int(count)
    specimens = list(specimen_runs)
    sequences = ["" for _ in sequence_index]
    for sequence, index in sequence_index.items():
        sequences[index] = sequence
    rows: list[int] = []
    columns: list[int] = []
    data: list[int] = []
    for row, specimen in enumerate(specimens):
        for column, count in specimen_rows[specimen].items():
            if count:
                rows.append(row)
                columns.append(column)
                data.append(count)
    counts = sparse.coo_matrix(
        (data, (rows, columns)), shape=(len(specimens), len(sequences)), dtype=np.int64
    ).tocsr()
    observations = pd.DataFrame(
        [
            {
                "project_id": project_id,
                "deposited_specimen_id": specimen,
                "source_runs": json.dumps(sorted(specimen_runs[specimen])),
            }
            for specimen in specimens
        ]
    )
    artifact = build_sv_artifact(
        counts,
        observations,
        sequences,
        dataset_id=dataset_id,
        tool_version="maliampi-tools",
        provenance={"source_rds": [str(path) for path in rds_paths]},
    )
    write_sv_h5ad(artifact, output_h5ad)
    # DADA2 expects sample-by-sequence integer data and exact dimnames; this is transient.
    write_rds(
        MatrixWrapper(counts.toarray().astype(np.int32), dimnames=(specimens, sequences)),
        str(output_rds),
    )
