# MaliAmPi

MaliAmPi is a modular, phylogenetic-placement amplicon workflow. Its canonical
sequence-variant contract is a versioned sparse H5AD artifact: observations are
project/specimen identities, variables are full sequence identities, and `X` is
an integer specimen-by-SV CSR count matrix.

The usual composed workflow is:

1. infer SVs with DADA2;
2. validate or build a reference package;
3. place with EPA-ng (default) or pplacer;
4. derive taxonomy and phylotypes in parallel;
5. calculate placement statistics.

Every stage is also independently runnable. Running raw FASTQ through every
terminal output is not required.

## Documentation

- [Architecture and H5AD contract](docs/architecture.md)
- [Manifest, identity, refpkg, and legacy-import inputs](docs/inputs.md)
- [Composed and standalone workflows](docs/workflows.md)
- [Resources, containers, CI, and releases](docs/operations.md)

## Requirements

- Nextflow 26.04.6 or newer
- Java 17 or newer
- Docker, or a configured Singularity environment

The `maliampi-tools` helper image is currently `linux/amd64` because its
`rds2py` dependency is not supported on ARM. Docker Desktop can run it under
emulation on Apple Silicon. CPU and memory are assigned through Nextflow process
labels in configuration, not hard-coded in workflow processes.

## Canonical SV artifact

`sv.h5ad` contains nonnegative integral counts in sparse CSR form, deposited
specimen and project identities in `obs`, normalized uppercase IUPAC sequences
and full SHA-256 identities in `var`, and schema/tool/refpkg provenance under
`uns["maliampi"]`. The helper package validates these invariants whenever an SV
artifact is read.

Legacy FASTA, long, map, weights, and share files are boundary imports or
optional terminal exports. They are not internal abundance contracts.

## Supported entrypoints

The composed workflow starts at `main.nf`. It accepts a manifest and can either
consume an existing refpkg or construct one. EPA-ng is the default placement
engine; pplacer is interchangeable for placement.

```sh
nextflow run main.nf -profile docker \
  --manifest manifest.csv \
  --project_id PROJECT --dataset_id DATASET \
  --refpkg refpkg.tar.gz --output output
```

Supported modular entrypoints are:

| Stage | Entrypoint |
|---|---|
| DADA2 SV inference | `modules/dada2.nf` |
| Swarm SV inference | `modules/swarm.nf` |
| refpkg construction | `modules/refpackage.nf` |
| placement | `subworkflows/local/place.nf -entry place_entry` |
| taxonomy | `subworkflows/local/taxonomy.nf -entry taxonomy_entry` |
| phylotypes | `subworkflows/local/phylotypes.nf -entry phylotypes_entry` |
| stats | `subworkflows/local/stats.nf -entry stats_entry` |
| legacy import/conversion | `subworkflows/local/convert.nf -entry normalize_sv_entry` |
| Good's filtering | `subworkflows/local/goods.nf` |

Taxonomy mappings and SV-to-phylotype mappings are Parquet. Rank/threshold
abundances are sparse H5AD with integer counts and a sparse relative-abundance
layer. Phylotype thresholds default to `0.1`, `0.3`, and `1.0`; all SVs must be
assigned, and new phylotypes may be added when required.

### Standalone Good's filtering

Good's filtering is opt-in and is never invoked by `main.nf`:

```sh
nextflow run subworkflows/local/goods.nf -profile docker \
  --sv_h5ad sv.h5ad --output output
```

Defaults are convergence delta `0.0001`, minimum reads `30`, minimum prevalence
`2`, seed `1`, and dropping non-converged specimens. Use
`--goods_keep_nonconverged true` to retain non-converged specimens; they never
contribute to prevalence. Outputs are `sv.goods.h5ad`,
`goods.convergence.parquet`, and `goods.curves.parquet`.

### Standalone Swarm

Swarm is an independently composable alternative SV-inference workflow and is
never invoked by `main.nf`:

```sh
nextflow run modules/swarm.nf -profile docker \
  --manifest manifest.csv \
  --project_id PROJECT --dataset_id DATASET --output output
```

vsearch performs read merging, filtering, dereplication, and chimera detection;
Swarm performs clustering. MaliAmPi owns registry parsing, sequence identities,
sparse count accumulation, and the final canonical `sv.h5ad`. The workflow also
publishes `swarm.stats.parquet` and the common specimen-failure report.

## Helper package development

The definitive internal helper code is the root `maliampi-tools` uv project.
Barcodecop is installed from PyPI; phylotypes remains an independent image and
release lifecycle.

```sh
uv sync --extra dev
uv lock --check
uv run ruff format --check src tests/python
uv run ruff check src tests/python
uv run ty check src
uv run pytest
```

Docker builds, Nextflow runs, image publication, pushes, and external archival
are explicit approval gates. The operational guide records the verified smoke
path and commands for reproducing it.

### Helper-image releases

The container workflow validates an AMD64 `maliampi-tools` build on relevant
pull requests and `master` changes without publishing. Pushing a `vX.Y.Z` tag
runs the same validation and then publishes the image to GHCR — version,
commit-SHA, and floating `latest` tags — and cuts a matching GitHub Release. See
[the container policy](docker/README.md).
