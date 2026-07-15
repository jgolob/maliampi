# Phase 06 — Committed config + centralized defaults

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/06-config` ·
**Files:** `nextflow.config`, `conf/base.config`, `conf/test.config`; CI version pin
only if required · **Depends on:** 01

## Goal

Create one authoritative configuration layer for resources, containers, and
defaults. Modules and imported workflows must not redeclare shared `params.*`.
Pin every documented and CI invocation to Nextflow **26.04.6**.

## Steps

1. Move existing `withLabel` resources into `conf/base.config`; preserve the
   current standard resource values and labels.
2. Put reduced resources and self-contained fixture inputs in `conf/test.config`.
   Use:
   ```nextflow
   params {
       manifest   = "${projectDir}/maliampi-practice-data/autotest_manifest.csv"
       repo_fasta = "${projectDir}/maliampi-practice-data/test.repo.fasta"
       repo_si    = "${projectDir}/maliampi-practice-data/test.repo.seq_info.csv"
       taxdmp     = "${projectDir}/maliampi-practice-data/taxdmp.zip"
       output     = 'output'
   }
   ```
3. Create profiles for `standard`, `docker`, `singularity`, `test`, and the
   backwards-compatible `testing` alias.
4. Declare `manifest.nextflowVersion = '>=26.04.6'`. Pin CI/setup examples to
   exactly `26.04.6`; do not retain the conflicting 24.10.0/26.04.3 examples.
5. Move shared scientific/container defaults out of module top levels as modules
   are touched in later phases. Phase 06 establishes the canonical values; do not
   duplicate them in newly split files.
6. Keep `nextflow.config.sample` temporarily; Phase 28 deprecates it in docs.

## Verify

```bash
NXF_VER=26.04.6 nextflow config -profile test,docker
NXF_VER=26.04.6 nextflow lint -o concise .
```

The resolved test manifest must exist, and lint output must not contain duplicate
parameter-definition warnings introduced by the new configuration.

## Commit / PR

- Title: `feat(config): committed config and Nextflow 26.04.6 baseline`
