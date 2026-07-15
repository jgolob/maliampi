# H5AD migration acceptance checklist

No Docker build or Nextflow smoke run should begin until the static gates below
are complete. CSV files may remain as terminal exports or immediate imports at
external-tool boundaries; they must not remain internal SV abundance contracts.

## Package contract

- [x] Versioned specimen-by-SV integer CSR H5AD creation and read validation.
- [x] Full sequence hashes and deterministic composite SV IDs.
- [x] Pure-Python RDS matrix serialization with integer storage and dimnames.
- [x] Official DADA2 `readRDS()` and `removeBimeraDenovo()` acceptance test.
- [x] H5AD-derived FASTA, gappa multiplicity, and optional pplacer exports.
- [x] Source H5AD file digest in derived abundance provenance.
- [ ] Legacy CSV/share/map/weights import that immediately emits H5AD.
- [ ] Refpkg SV registry creation, validation, and extension.
- [ ] Collision, duplicated-run, empty/rejected, large-count, and malformed-input tests.
- [ ] Scalability benchmark against the legacy sparse combiner.

## Active workflows

- [ ] Barcodecop runs from the unified helper image installed from PyPI.
- [ ] DADA2 emits only transient RDS plus pre/post-chimera H5AD internally.
- [ ] Placement accepts H5AD and creates only boundary inputs.
- [ ] Taxonomy mapping is Parquet and rank abundance artifacts are sparse H5AD.
- [ ] Phylotype mapping is Parquet and threshold abundance artifacts are sparse H5AD.
- [ ] Conversion accepts legacy formats only as import formats.
- [ ] Refpkg validates full SV identities before placement.
- [ ] Stats and refpkg construction no longer depend on retired helper images.

## Legacy and cleanup

- [ ] Swarm uses H5AD internally.
- [ ] Both legacy placement monoliths use H5AD internally or are explicitly retired.
- [ ] Good's filter process, parameters, image, and code are removed.
- [ ] Retired DADA2 conversion processes and custom image parameters are removed.
- [ ] `container__fastatools`, `container__dada2pplacer`, and other superseded helpers are gone.
- [ ] `vendor/` and duplicate imported Docker directories are removed only after replacement review.

## Verification sequence

- [ ] Ruff format/lint, ty, pytest, and `uv lock --check` pass.
- [ ] Nextflow configuration and scripts parse without undefined-container warnings.
- [ ] Golden legacy exports match H5AD-derived exports.
- [ ] AMD64 helper image builds once and its container-level acceptance tests pass.
- [ ] Sample workflow succeeds from FASTQ through canonical `sv.h5ad`.
- [ ] Full workflow succeeds through placement, taxonomy, phylotypes, and terminal outputs.
