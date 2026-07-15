# Phase 25 — Core per-entry CI and monolith-job migration

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/25-entry-ci` ·
**Files:** `.github/workflows/test.yaml`, assertion helpers/fixtures as needed ·
**Depends on:** 09a, 12a, 13–18

## Goal

Prove the core modular entries are supported independently and migrate every CI job
away from the old monoliths before Phase 19 deletes them. Pin Nextflow exactly to
**26.04.6** in every job.

## Required jobs

1. `validate_refpkg_entry`: valid fixture plus tampered-package rejection.
2. Conversion named entries: SV normalization modes and alignment preparation.
3. `sv_entry`: canonical four-file output.
4. `build_refpkg_entry`: archive, validation report, and stable `refpkg_id`.
5. `place_entry` matrix: EPA-ng and pplacer; deduplicated plus specimen jplace.
6. `taxonomy_entry`: per-SV and count/relative-abundance tables.
7. `stats_entry`: modern default; a smaller explicit legacy case.
8. Full composition with provided refpkg: build stage absent, peer branches present.
9. Skip-flag matrix for the branches available before Phase 22.

Use Phase-24 helpers for content, identity, and count-conservation assertions.

## Coordination

Remove or rewrite all jobs that invoke `epang_place_classify.nf` or
`pplacer_place_classify.nf`. Phase 19 depends on this phase being merged and green.
Phylotypes/set-extension jobs are added later in Phase 25a so they do not block the
safe deletion ordering.

## Verify

- All jobs use Nextflow 26.04.6.
- Repository search finds no live CI reference to either monolith.
- Corrupting a canonical output or `refpkg_id` fails CI.

## Commit / PR

- Title: `test(ci): core modular entries and monolith-job migration`
