"""Pipeline-facing command-line adapters for maliampi-tools."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .aggregate import (
    aggregate_mapping,
    aggregate_taxonomy_rank,
    write_phylotype_mapping,
    write_taxonomy_mapping,
)
from .combine import combine_seqtabs
from .export import export_legacy
from .goods import write_goods_outputs
from .importers import import_legacy_sv
from .refpkg import validate_sv_registry, write_sv_registry
from .sv import filter_sv_artifact, read_sv_h5ad, sequence_hash, write_sv_h5ad
from .swarm import (
    combine_dsv_registries,
    filter_swarm_seeds,
    finalize_swarm_h5ad,
    prepare_specimen_dsv,
)


def combine_seqtabs_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--output-rds", required=True)
    parser.add_argument("rds", nargs="+")
    args = parser.parse_args()
    combine_seqtabs(
        args.rds,
        args.output_h5ad,
        args.output_rds,
        project_id=args.project_id,
        dataset_id=args.dataset_id,
    )


def filter_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_h5ad")
    parser.add_argument("retained_hashes")
    parser.add_argument("output_h5ad")
    parser.add_argument("--sequences", action="store_true")
    args = parser.parse_args()
    retained = [
        line.strip() for line in Path(args.retained_hashes).read_text().splitlines() if line.strip()
    ]
    if args.sequences:
        retained = [sequence_hash(sequence) for sequence in retained]
    artifact = read_sv_h5ad(args.input_h5ad)
    write_sv_h5ad(filter_sv_artifact(artifact, retained), args.output_h5ad)


def validate_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad")
    args = parser.parse_args()
    read_sv_h5ad(args.h5ad)


def export_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad")
    parser.add_argument("output_dir")
    parser.add_argument("--pplacer-reduplication", action="store_true")
    args = parser.parse_args()
    export_legacy(args.h5ad, args.output_dir, pplacer_reduplication=args.pplacer_reduplication)


def aggregate_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sv_h5ad")
    parser.add_argument("mapping_parquet")
    parser.add_argument("output_h5ad")
    parser.add_argument("--identity-column", required=True)
    parser.add_argument("--metadata-json", default="{}")
    parser.add_argument("--placement-identity")
    parser.add_argument("--refpkg-identity")
    args = parser.parse_args()
    metadata = json.loads(args.metadata_json)
    if args.placement_identity:
        metadata["placement_identity"] = json.loads(Path(args.placement_identity).read_text())
    if args.refpkg_identity:
        metadata["refpkg_identity"] = json.loads(Path(args.refpkg_identity).read_text())
    aggregate_mapping(
        args.sv_h5ad,
        args.mapping_parquet,
        args.output_h5ad,
        identity_column=args.identity_column,
        metadata=metadata,
    )


def phylotype_mapping_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("parquet")
    args = parser.parse_args()
    write_phylotype_mapping(args.csv, args.parquet)


def taxonomy_mapping_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("parquet")
    parser.add_argument("--sv-h5ad")
    args = parser.parse_args()
    write_taxonomy_mapping(args.csv, args.parquet, sv_h5ad=args.sv_h5ad)


def taxonomy_aggregate_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sv_h5ad")
    parser.add_argument("taxonomy_parquet")
    parser.add_argument("output_h5ad")
    parser.add_argument("--rank", required=True)
    parser.add_argument("--placement-identity")
    parser.add_argument("--refpkg-identity")
    args = parser.parse_args()
    metadata: dict[str, object] = {}
    if args.placement_identity:
        metadata["placement_identity"] = json.loads(Path(args.placement_identity).read_text())
    if args.refpkg_identity:
        metadata["refpkg_identity"] = json.loads(Path(args.refpkg_identity).read_text())
    aggregate_taxonomy_rank(
        args.sv_h5ad,
        args.taxonomy_parquet,
        args.output_h5ad,
        rank=args.rank,
        metadata=metadata,
    )


def import_legacy_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--long")
    parser.add_argument("--sharetable")
    parser.add_argument("--sv-map")
    parser.add_argument("--sv-weights")
    args = parser.parse_args()
    import_legacy_sv(
        args.fasta,
        args.output_h5ad,
        project_id=args.project_id,
        dataset_id=args.dataset_id,
        long=args.long,
        sharetable=args.sharetable,
        sv_map=args.sv_map,
        sv_weights=args.sv_weights,
    )


def registry_write_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sv_h5ad")
    parser.add_argument("registry_parquet")
    args = parser.parse_args()
    write_sv_registry(args.sv_h5ad, args.registry_parquet)


def registry_validate_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sv_h5ad")
    parser.add_argument("registry_parquet")
    args = parser.parse_args()
    print(json.dumps(validate_sv_registry(args.sv_h5ad, args.registry_parquet), sort_keys=True))


def goods_filter_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_h5ad")
    parser.add_argument("output_h5ad")
    parser.add_argument("--convergence-parquet", default="goods.convergence.parquet")
    parser.add_argument("--curves-parquet", default="goods.curves.parquet")
    parser.add_argument("--convergence-delta", type=float, default=0.0001)
    parser.add_argument("--min-reads", type=int, default=30)
    parser.add_argument("--min-prevalence", type=int, default=2)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--keep-nonconverged", action="store_true")
    args = parser.parse_args()
    write_goods_outputs(
        args.input_h5ad,
        args.output_h5ad,
        args.convergence_parquet,
        args.curves_parquet,
        convergence_delta=args.convergence_delta,
        min_reads=args.min_reads,
        min_prevalence=args.min_prevalence,
        seed=args.seed,
        keep_nonconverged=args.keep_nonconverged,
    )


def swarm_specimen_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("derep_fasta")
    parser.add_argument("output_registry")
    parser.add_argument("--specimen", required=True)
    parser.add_argument("--batch", required=True)
    args = parser.parse_args()
    prepare_specimen_dsv(
        args.derep_fasta,
        args.output_registry,
        specimen=args.specimen,
        batch=args.batch,
    )


def swarm_combine_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-registry", required=True)
    parser.add_argument("registries", nargs="+")
    args = parser.parse_args()
    combine_dsv_registries(args.registries, args.output_fasta, args.output_registry)


def swarm_filter_seeds_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("output_fasta")
    args = parser.parse_args()
    filter_swarm_seeds(args.input_fasta, args.output_fasta)


def swarm_finalize_main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("registry_parquet")
    parser.add_argument("clusters")
    parser.add_argument("nonchimera_fasta")
    parser.add_argument("output_h5ad")
    parser.add_argument("stats_parquet")
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--dataset-id", required=True)
    args = parser.parse_args()
    finalize_swarm_h5ad(
        args.registry_parquet,
        args.clusters,
        args.nonchimera_fasta,
        args.output_h5ad,
        args.stats_parquet,
        project_id=args.project_id,
        dataset_id=args.dataset_id,
    )
