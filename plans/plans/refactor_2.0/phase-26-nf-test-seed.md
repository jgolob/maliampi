# Phase 26 — Seed nf-test

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/26-nf-test` ·
**Files:** `nf-test.config`, `tests/**/*.nf.test`, CI ·
**Depends on:** 25a and the entries under test

## Goal

Introduce nf-test around functions and modular contracts while retaining the small
end-to-end CI jobs. Pin its Nextflow runtime to 26.04.6.

## Tests

- Test `read_manifest` as a function or through a tiny test wrapper workflow; it is
  not a process and must not use a `nextflow_process` skeleton.
- Test canonical SV normalization for each accepted source mode and invalid mixes.
- Test refpkg validation identity and representative failures.
- Test placement subworkflow outputs for EPA-ng; use a separate pplacer fixture test
  where runtime permits.
- Test taxonomy and phylotype abundance conservation.
- Test effective threshold parsing (replace, append, sort, deduplicate, invalid).
- Test phylotype-set extension invariants: stable old IDs, join-old, create-new,
  add-threshold, complete mappings, new set ID, and parent linkage.

Prefer structural/content assertions over large opaque snapshots. Snapshot small,
deterministic CSV/JSON only; jplace uses semantic assertions.

## Verify

```bash
nf-test test
```

Mutate a canonical ID, checksum, threshold mapping, and count total in turn and
confirm each relevant test fails. Never update snapshots blindly.

## Commit / PR

- Title: `test: nf-test modular contracts and versioned-set extension`
