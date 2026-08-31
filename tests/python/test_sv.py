from __future__ import annotations

from hashlib import sha256

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from maliampi_tools.aggregate import (
    aggregate_mapping,
    aggregate_taxonomy_rank,
    write_phylotype_mapping,
    write_taxonomy_mapping,
)
from maliampi_tools.export import export_legacy
from maliampi_tools.importers import import_legacy_sv
from maliampi_tools.refpkg import validate_sv_registry, write_sv_registry
from maliampi_tools.sv import SvValidationError, build_sv_artifact, read_sv_h5ad, write_sv_h5ad


def artifact():
    return build_sv_artifact(
        sparse.csr_matrix([[2, 0], [3, 5]], dtype=np.int64),
        pd.DataFrame(
            {
                "project_id": ["project-a", "project-a"],
                "deposited_specimen_id": ["sample-1", "sample-2"],
                "batches": ["run-1", "run-2"],
            }
        ),
        ["acgt", "ACGTR"],
        dataset_id="dataset-a",
        tool_version="test",
    )


def test_h5ad_contract_round_trip_and_exports(tmp_path):
    path = tmp_path / "sv.h5ad"
    write_sv_h5ad(artifact(), path)
    restored = read_sv_h5ad(path)
    assert sparse.isspmatrix_csr(restored.X)
    assert restored.X.dtype == np.int64
    assert restored.var["sequence"].tolist() == ["ACGT", "ACGTR"]
    export_legacy(path, tmp_path / "exports", pplacer_reduplication=True)
    assert (tmp_path / "exports" / "sv.fasta").read_text().count(">sv-") == 2
    assert (tmp_path / "exports" / "sv.map.csv").exists()


def test_bad_sequence_and_orphan_are_rejected():
    with pytest.raises(SvValidationError, match="Invalid IUPAC"):
        build_sv_artifact(
            [[1]],
            pd.DataFrame({"project_id": ["p"], "deposited_specimen_id": ["s"]}),
            ["ACGTZ"],
            dataset_id="d",
        )
    with pytest.raises(SvValidationError, match="orphan"):
        build_sv_artifact(
            [[0]],
            pd.DataFrame({"project_id": ["p"], "deposited_specimen_id": ["s"]}),
            ["ACGT"],
            dataset_id="d",
        )


def test_complete_parquet_mapping_preserves_counts_and_relative_abundance(tmp_path):
    source = tmp_path / "sv.h5ad"
    write_sv_h5ad(artifact(), source)
    mapping = pd.DataFrame(
        {
            "sv_id": artifact().var_names,
            "phylotype": ["pt-1", "pt-1"],
        }
    )
    mapping_path = tmp_path / "mapping.parquet"
    mapping.to_parquet(mapping_path)
    output = tmp_path / "phylotypes.h5ad"
    aggregate_mapping(
        source, mapping_path, output, identity_column="phylotype", metadata={"threshold": 0.1}
    )
    result = ad.read_h5ad(output)
    assert result.X.sum() == 10
    assert np.allclose(np.asarray(result.layers["relative_abundance"].sum(axis=1)).ravel(), 1.0)
    assert (
        result.uns["maliampi"]["source_sv_artifact_digest"]
        == sha256(source.read_bytes()).hexdigest()
    )


def test_external_mappings_are_imported_to_parquet_and_aggregate_from_h5ad(tmp_path):
    source = tmp_path / "sv.h5ad"
    write_sv_h5ad(artifact(), source)
    phylotype_csv = tmp_path / "phylotypes.csv"
    pd.DataFrame({"phylotype": ["pt-1", "pt-2"], "sv": artifact().var_names}).to_csv(
        phylotype_csv, index=False
    )
    phylotype_parquet = tmp_path / "phylotypes.parquet"
    write_phylotype_mapping(phylotype_csv, phylotype_parquet)
    assert pd.read_parquet(phylotype_parquet)["sv_id"].tolist() == artifact().var_names.tolist()

    taxonomy_csv = tmp_path / "taxonomy.csv"
    pd.DataFrame(
        {
            "sv": artifact().var_names,
            "want_rank": ["genus", "genus"],
            "rank": ["genus", "genus"],
            "tax_name": ["Genus A", None],
        }
    ).to_csv(taxonomy_csv, index=False)
    taxonomy_parquet = tmp_path / "taxonomy.parquet"
    write_taxonomy_mapping(taxonomy_csv, taxonomy_parquet, sv_h5ad=source)
    output = tmp_path / "genus.h5ad"
    provenance = {
        "placement_identity": {"placer": "epang", "refpkg_id": "refpkg-1"},
        "refpkg_identity": {"refpkg_id": "refpkg-1"},
    }
    aggregate_taxonomy_rank(source, taxonomy_parquet, output, rank="genus", metadata=provenance)
    result = ad.read_h5ad(output)
    assert set(result.var_names) == {"Genus A", "unclassified"}
    assert result.X.sum() == 10
    assert "assigned_rank" in pd.read_parquet(taxonomy_parquet).columns
    assert result.uns["maliampi"]["placement_identity"]["placer"] == "epang"
    assert result.uns["maliampi"]["refpkg_identity"]["refpkg_id"] == "refpkg-1"

    distant_csv = tmp_path / "distant.csv"
    pd.DataFrame(
        {"sv": [artifact().var_names[0]], "want_rank": ["species"], "tax_name": ["Species A"]}
    ).to_csv(distant_csv, index=False)
    distant_parquet = tmp_path / "distant.parquet"
    write_taxonomy_mapping(distant_csv, distant_parquet, sv_h5ad=source)
    imported = pd.read_parquet(distant_parquet)
    assert set(imported.loc[imported["rank"] == "species", "sv_id"]) == set(artifact().var_names)


def test_legacy_long_import_immediately_creates_a_valid_h5ad(tmp_path):
    fasta = tmp_path / "legacy.fasta"
    fasta.write_text(">old-1\nACGT\n>old-2\nACGTR\n")
    long = tmp_path / "legacy.long.csv"
    pd.DataFrame(
        {
            "specimen": ["sample-1", "sample-1", "sample-2"],
            "sv": ["old-1", "old-1", "old-2"],
            "count": [2, 3, 5],
        }
    ).to_csv(long, index=False)
    output = tmp_path / "sv.h5ad"
    import_legacy_sv(
        fasta,
        output,
        project_id="project-a",
        dataset_id="dataset-a",
        long=long,
    )
    imported = read_sv_h5ad(output)
    assert imported.X.sum() == 10
    assert imported.n_obs == 2
    assert imported.n_vars == 2


def test_legacy_long_import_gzipped(tmp_path):
    """Gzipped FASTA and long CSV are transparently decompressed."""
    import gzip

    fasta_gz = tmp_path / "legacy.fasta.gz"
    with gzip.open(fasta_gz, "wt") as f:
        f.write(">old-1\nACGT\n>old-2\nACGTR\n")
    long_gz = tmp_path / "legacy.long.csv.gz"
    with gzip.open(long_gz, "wt") as f:
        f.write("specimen,sv,count\nsample-1,old-1,2\nsample-1,old-1,3\nsample-2,old-2,5\n")
    output = tmp_path / "sv.h5ad"
    import_legacy_sv(
        fasta_gz,
        output,
        project_id="project-a",
        dataset_id="dataset-a",
        long=long_gz,
    )
    imported = read_sv_h5ad(output)
    assert imported.X.sum() == 10
    assert imported.n_obs == 2
    assert imported.n_vars == 2


def test_refpkg_registry_preserves_full_hash_and_composite_id_identity(tmp_path):
    source = tmp_path / "sv.h5ad"
    write_sv_h5ad(artifact(), source)
    registry = tmp_path / "sv-registry.parquet"
    write_sv_registry(source, registry)
    identity = validate_sv_registry(source, registry)
    assert len(identity["registry_sha256"]) == 64

    broken = pd.read_parquet(registry)
    broken.loc[0, "sv_id"] = "sv-0000000000-0000000000000000"
    broken.to_parquet(registry, index=False)
    with pytest.raises(SvValidationError, match="composite SV ID"):
        validate_sv_registry(source, registry)
