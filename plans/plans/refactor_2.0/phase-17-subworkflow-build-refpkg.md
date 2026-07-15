# Phase 17 — Standalone build-and-validate-refpkg subworkflow

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/17-subwf-build-refpkg` ·
**Files:** `subworkflows/local/build_refpkg.nf`, `modules/refpackage.nf` signature ·
**Depends on:** 09a

## Goal

Build a refpkg from explicit inputs, validate it immediately, and emit the immutable
archive together with its identity.

## Contract

```text
build_refpkg(sv_fasta, repo_fasta, repo_seq_info,
             taxdmp=null, covariance_model=null)
```

The imported workflow must not read repository inputs from `params`. The standalone
entry translates `--sv_fasta`, `--repo_fasta`, `--repo_si`, optional `--taxdmp`, and
optional `--rfam` into channels.

Emit `refpkg/refpkg.tar.gz`, `refpkg/refpkg.identity.json`, and the validation
report. Do not emit an unvalidated package.

## Verify

- Fixture build produces a valid refpkg and stable `refpkg_id`.
- Repacking identical logical content preserves the ID.
- Taxonomy download credentials/email are required only if that branch actually
  performs a download.

## Commit / PR

- Title: `feat(refpkg): explicit build inputs and validated identity output`
