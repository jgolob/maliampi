"""Sparse taxonomy and phylotype abundance aggregation."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .sv import SvValidationError, read_sv_h5ad

TAXONOMY_RANKS = ("superkingdom", "phylum", "class", "order", "family", "genus", "species")


def write_phylotype_mapping(csv_path: str | Path, output_parquet: str | Path) -> None:
    """Import the external phylotypes CSV into the internal Parquet contract."""
    mapping = pd.read_csv(csv_path)
    if set(mapping.columns) != {"phylotype", "sv"}:
        raise SvValidationError("Phylotype mapping must have exactly phylotype and sv columns")
    mapping = mapping.rename(columns={"sv": "sv_id"})
    if mapping["sv_id"].isna().any() or mapping["phylotype"].isna().any():
        raise SvValidationError("Phylotype mapping contains an unassigned SV")
    if mapping["sv_id"].duplicated().any():
        raise SvValidationError("Phylotype mapping assigns an SV more than once")
    mapping.to_parquet(output_parquet, index=False)


def write_taxonomy_mapping(
    csv_path: str | Path,
    output_parquet: str | Path,
    *,
    sv_h5ad: str | Path | None = None,
) -> None:
    """Import gappa taxonomy output into a complete rank-specific Parquet mapping."""
    mapping = pd.read_csv(csv_path)
    required = {"sv", "want_rank", "tax_name"}
    if not required.issubset(mapping.columns):
        raise SvValidationError(f"Taxonomy mapping must contain {sorted(required)}")
    # Gappa's normalized table retains its assigned rank as ``rank`` while
    # ``want_rank`` names the abundance level to construct.  Preserve both
    # rather than creating duplicate dataframe labels during the rename.
    if "rank" in mapping.columns:
        mapping = mapping.rename(columns={"rank": "assigned_rank"})
    mapping = mapping.rename(columns={"sv": "sv_id", "want_rank": "rank"})
    mapping["taxon"] = mapping["tax_name"].fillna("unclassified").astype(str)
    if mapping.duplicated(["sv_id", "rank"]).any():
        raise SvValidationError("Taxonomy mapping assigns an SV more than once at a rank")
    if sv_h5ad is not None:
        source = read_sv_h5ad(sv_h5ad)
        ranks = tuple(sorted(set(mapping["rank"]))) or TAXONOMY_RANKS
        missing = [
            {"sv_id": sv_id, "rank": rank, "taxon": "unclassified"}
            for rank in ranks
            for sv_id in source.var_names
            if not ((mapping["rank"] == rank) & (mapping["sv_id"] == sv_id)).any()
        ]
        if missing:
            mapping = pd.concat((mapping, pd.DataFrame(missing)), ignore_index=True)
    mapping.to_parquet(output_parquet, index=False)


def aggregate_taxonomy_rank(
    sv_h5ad: str | Path,
    taxonomy_parquet: str | Path,
    output_h5ad: str | Path,
    *,
    rank: str,
    metadata: dict[str, object] | None = None,
) -> None:
    """Create a count and sparse relative-abundance artifact for one taxonomy rank."""
    mapping = pd.read_parquet(taxonomy_parquet)
    selected = mapping.loc[mapping["rank"] == rank, ["sv_id", "taxon"]]
    temporary_mapping = Path(output_h5ad).with_suffix(".mapping.parquet")
    selected.to_parquet(temporary_mapping, index=False)
    try:
        aggregate_mapping(
            sv_h5ad,
            temporary_mapping,
            output_h5ad,
            identity_column="taxon",
            metadata={"rank": rank, **(metadata or {})},
        )
    finally:
        temporary_mapping.unlink(missing_ok=True)


def aggregate_mapping(
    sv_h5ad: str | Path,
    mapping_parquet: str | Path,
    output_h5ad: str | Path,
    *,
    identity_column: str,
    metadata: dict[str, object],
) -> None:
    """Aggregate an SV matrix by a complete one-to-one Parquet mapping."""
    source_path = Path(sv_h5ad)
    source = read_sv_h5ad(source_path)
    mapping = pd.read_parquet(mapping_parquet)
    if {"sv_id", identity_column} - set(mapping.columns):
        raise SvValidationError(f"Mapping must contain sv_id and {identity_column}")
    if mapping["sv_id"].duplicated().any() or set(mapping["sv_id"]) != set(source.var_names):
        raise SvValidationError("Mapping must assign every SV exactly once")
    target_ids = pd.Index(pd.unique(mapping[identity_column]), name=identity_column)
    target_index = {identifier: index for index, identifier in enumerate(target_ids)}
    columns = mapping.set_index("sv_id").loc[source.var_names, identity_column]
    transform = sparse.csr_matrix(
        (
            np.ones(source.n_vars, dtype=np.int8),
            (np.arange(source.n_vars), [target_index[value] for value in columns]),
        ),
        shape=(source.n_vars, len(target_ids)),
    )
    counts = (source.X @ transform).tocsr().astype(np.int64)
    row_totals = np.asarray(counts.sum(axis=1)).ravel()
    relative = (
        sparse.diags(
            np.divide(
                1.0, row_totals, out=np.zeros_like(row_totals, dtype=float), where=row_totals != 0
            )
        )
        @ counts
    )
    var = pd.DataFrame({identity_column: target_ids}, index=target_ids)
    # Rank and threshold are per-feature scalar annotations.  Provenance is
    # structured and belongs only in ``uns["maliampi"]``; writing a dict as a
    # dataframe column is invalid for H5AD serialization.
    for key, value in metadata.items():
        if value is None or isinstance(value, str | int | float | bool):
            var[key] = value
    result = ad.AnnData(X=counts, obs=cast(pd.DataFrame, source.obs).copy(), var=var)
    result.layers["relative_abundance"] = relative.tocsr()
    result.uns["maliampi"] = {
        "source_sv_artifact_digest": sha256(source_path.read_bytes()).hexdigest(),
        "source_dataset_id": source.uns["maliampi"]["dataset_id"],
        **metadata,
    }
    cast(Any, result).write_h5ad(str(output_h5ad))
