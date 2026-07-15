# Phase 14 — Standalone taxonomy subworkflow

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/14-subwf-taxonomy` ·
**Files:** `subworkflows/local/taxonomy.nf`, `modules/taxonomy.nf` publication paths ·
**Depends on:** 09a, 11

## Goal

Create the common gappa taxonomy peer pathway. It consumes placement from either
engine and remains entirely independent of phylotypes.

## Contract

```text
taxonomy(dedup_jplace, validated_refpkg, sv_long)
```

Standalone `taxonomy_entry` accepts `--dedup_jplace`, `--refpkg`, and `--sv_long`.
It emits under `taxonomy/`:

- Per-SV taxonomy in long form
- Per-rank specimen count tables
- Per-rank specimen relative-abundance tables

## Steps

1. Extract taxonomy/refpkg components from the validated package.
2. Run gappa assignment and extract the per-SV taxonomy.
3. Join the per-SV assignments to `sv_long` for each supported rank.
4. Change old `classify/` publication paths to the intentional clean v2
   `taxonomy/` layout.
5. Do not include or call phylotypes. Do not branch on the placer.

## Verify

- EPA-ng and pplacer fixture placements both traverse the same taxonomy workflow.
- Per-SV taxonomy and every wide table are nonempty.
- Per-specimen count totals are conserved; nonempty relative-abundance rows sum to
  1 within numeric tolerance.

## Commit / PR

- Title: `feat(taxonomy): standalone common gappa taxonomy pathway`
