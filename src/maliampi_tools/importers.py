"""One-way import of legacy SV abundance formats into canonical H5AD."""

from __future__ import annotations

import csv
import gzip
from collections import defaultdict
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import IO

import numpy as np
import pandas as pd
from scipy import sparse

from .sv import SvValidationError, build_sv_artifact, write_sv_h5ad


@contextmanager
def _open_maybe_gzipped(path: str | Path, *, newline: str | None = "") -> Iterator[IO[str]]:
    """Open a file for reading text, transparently decompressing .gz files."""
    p = Path(path)
    handle = gzip.open(p, "rt", newline=newline) if p.suffix == ".gz" else p.open(newline=newline)  # noqa: SIM115
    try:
        yield handle
    finally:
        handle.close()


def _fasta_records(path: str | Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence: list[str] = []
    with _open_maybe_gzipped(path, newline=None) as handle:
        for raw_line in handle.read().splitlines():
            if raw_line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence)))
                identifier = raw_line[1:].split(maxsplit=1)[0]
                sequence = []
            else:
                sequence.append(raw_line.strip())
        if identifier is not None:
            records.append((identifier, "".join(sequence)))
    if not records or len({identifier for identifier, _ in records}) != len(records):
        raise SvValidationError("Legacy FASTA must contain one or more uniquely named sequences")
    return records


def _long_rows(path: str | Path) -> list[tuple[str, str, int]]:
    rows: list[tuple[str, str, int]] = []
    with _open_maybe_gzipped(path) as handle:
        reader = csv.DictReader(handle)
        if not {"specimen", "sv", "count"}.issubset(reader.fieldnames or []):
            raise SvValidationError("Legacy long table requires specimen, sv, and count columns")
        for row in reader:
            try:
                count = int(row["count"])
            except (TypeError, ValueError) as error:
                raise SvValidationError("Legacy long table has a non-integer count") from error
            rows.append((row["specimen"], row["sv"], count))
    return rows


def _share_rows(path: str | Path) -> list[tuple[str, str, int]]:
    rows: list[tuple[str, str, int]] = []
    with _open_maybe_gzipped(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        group = next((field for field in fields if field.lower() == "group"), None)
        reserved = {"label", "group", "numotus", "numsvs"}
        variants = [field for field in fields if field.lower() not in reserved]
        if group is None or not variants:
            raise SvValidationError("Legacy share table requires group and SV columns")
        for row in reader:
            for variant in variants:
                try:
                    count = int(row[variant])
                except (TypeError, ValueError) as error:
                    raise SvValidationError("Legacy share table has a non-integer count") from error
                rows.append((row[group], variant, count))
    return rows


def _map_weight_rows(map_path: str | Path, weights_path: str | Path) -> list[tuple[str, str, int]]:
    specimen_for: dict[str, str] = {}
    with _open_maybe_gzipped(map_path) as handle:
        for row in csv.reader(handle):
            if len(row) != 2 or not all(row):
                raise SvValidationError("Legacy map must be headerless specimen-SV ID,specimen")
            specimen_for[row[0]] = row[1]
    rows: list[tuple[str, str, int]] = []
    with _open_maybe_gzipped(weights_path) as handle:
        for row in csv.reader(handle):
            if len(row) != 3 or row[1] not in specimen_for:
                raise SvValidationError("Legacy weights/map pair is inconsistent")
            try:
                count = int(row[2])
            except ValueError as error:
                raise SvValidationError("Legacy weights table has a non-integer count") from error
            rows.append((specimen_for[row[1]], row[0], count))
    return rows


def import_legacy_sv(
    fasta: str | Path,
    output_h5ad: str | Path,
    *,
    project_id: str,
    dataset_id: str,
    long: str | Path | None = None,
    sharetable: str | Path | None = None,
    sv_map: str | Path | None = None,
    sv_weights: str | Path | None = None,
) -> None:
    """Import exactly one legacy abundance representation into canonical H5AD."""
    modes = [long is not None, sharetable is not None, sv_map is not None or sv_weights is not None]
    if sum(modes) != 1 or ((sv_map is None) != (sv_weights is None)):
        raise SvValidationError("Provide exactly one of long, sharetable, or a map/weights pair")
    records = _fasta_records(fasta)
    sequence_for = dict(records)
    if long is not None:
        rows = _long_rows(long)
    elif sharetable is not None:
        rows = _share_rows(sharetable)
    else:
        if sv_map is None or sv_weights is None:  # guarded above; narrows static type.
            raise AssertionError("map/weights mode requires both inputs")
        rows = _map_weight_rows(sv_map, sv_weights)
    counts_by_sample: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for specimen, sv_id, count in rows:
        if not specimen or sv_id not in sequence_for or count < 0:
            raise SvValidationError("Legacy abundance data has an invalid specimen, SV, or count")
        counts_by_sample[specimen][sv_id] += count
    specimens = sorted(counts_by_sample)
    legacy_ids = [identifier for identifier, _ in records]
    columns = {identifier: index for index, identifier in enumerate(legacy_ids)}
    data: list[int] = []
    row_indices: list[int] = []
    column_indices: list[int] = []
    for row_index, specimen in enumerate(specimens):
        for sv_id, count in counts_by_sample[specimen].items():
            if count:
                row_indices.append(row_index)
                column_indices.append(columns[sv_id])
                data.append(count)
    matrix = sparse.coo_matrix(
        (data, (row_indices, column_indices)),
        shape=(len(specimens), len(records)),
        dtype=np.int64,
    ).tocsr()
    observations = pd.DataFrame(
        {
            "project_id": project_id,
            "deposited_specimen_id": specimens,
            "source_format": "legacy-import",
        }
    )
    artifact = build_sv_artifact(
        matrix,
        observations,
        [sequence for _, sequence in records],
        dataset_id=dataset_id,
        tool_version="maliampi-tools",
        provenance={"legacy_fasta": str(fasta)},
    )
    write_sv_h5ad(artifact, output_h5ad)
