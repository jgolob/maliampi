# Workflows

All commands below use Docker. Replace `-profile docker` with the configured
Singularity profile when appropriate, and supply a separate Nextflow work
directory for durable runs.

## Composed workflow

This is the normal raw-FASTQ path. EPA-ng is the default placement engine;
`--placer pplacer` selects the interchangeable alternative.

```sh
nextflow run main.nf -profile docker \
  --manifest manifest.csv \
  --project_id PROJECT --dataset_id DATASET \
  --refpkg existing-refpkg.tar.gz \
  --output output
```

Omit `--refpkg` and provide `--repo_fasta`, `--repo_si`, and `--email` to build
a new package. Use `--sv_only true` to stop after canonical SV generation.
`--skip_taxonomy`, `--skip_phylotypes`, and `--skip_stats` selectively omit
terminal analyses.

The composed run publishes `output/sv/sv.h5ad`, a refpkg when constructed,
placement files, taxonomy Parquet/H5AD artifacts, phylotype Parquet/H5AD
artifacts, and placement statistics.

## Standalone Good's filtering

Good's filtering is opt-in and never runs from `main.nf`.

```sh
nextflow run subworkflows/local/goods.nf -profile docker \
  --sv_h5ad sv.h5ad --output output
```

It publishes `sv.goods.h5ad`, `goods.convergence.parquet`, and
`goods.curves.parquet`. Defaults are delta `0.0001`, minimum reads `30`,
minimum prevalence `2`, and seed `1`. `--goods_keep_nonconverged true` retains
non-converged specimens, but they never contribute to SV prevalence.

## Standalone Swarm

Swarm is an independent FASTQ-to-H5AD pathway, not an engine option in
`main.nf`.

```sh
nextflow run modules/swarm.nf -profile docker \
  --manifest manifest.csv \
  --project_id PROJECT --dataset_id DATASET \
  --output output
```

vsearch handles merging, filtering, dereplication, and chimera removal; Swarm
does clustering. MaliAmPi creates the DSV registry, applies canonical
identities, accumulates sparse counts, and publishes `sv.h5ad` plus
`swarm.stats.parquet`.

## Modular downstream stages

| Stage | Entrypoint | Required durable inputs |
|---|---|---|
| Build refpkg | `modules/refpackage.nf` | `sv_fasta`, `sv_h5ad`, repository inputs |
| Placement | `subworkflows/local/place.nf -entry place_entry` | `sv_h5ad`, `refpkg` |
| Taxonomy | `subworkflows/local/taxonomy.nf -entry taxonomy_entry` | `dedup_jplace`, `refpkg`, `sv_h5ad`, placement identity |
| Phylotypes | `subworkflows/local/phylotypes.nf -entry phylotypes_entry` | `dedup_jplace`, `sv_h5ad`, placement identity, dataset ID |
| Stats | `subworkflows/local/stats.nf -entry stats_entry` | deduplicated jplace and specimen jplaces |
| Legacy import | `subworkflows/local/convert.nf -entry normalize_sv_entry` | one supported legacy input set |

Phylotypes default to thresholds `0.1`, `0.3`, and `1.0`. Add values with
`--phylotype_add_thresholds`; use calibrated
`--phylotype_kr_thresholds` only with `--phylotype_distance kr`.
