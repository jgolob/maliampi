# Phase 15 — Standalone modern/legacy stats subworkflow

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/15-subwf-stats` ·
**Files:** `subworkflows/local/stats.nf`, `modules/stats.nf` publication paths ·
**Depends on:** 12, 13

## Goal

Create the stats peer pathway with modern gappa metrics by default and preserved
legacy guppy/pplacer metrics only by explicit request.

## Contract

```text
stats(dedup_jplace, specimen_jplaces,
      legacy_stats=false, validated_refpkg=null, sv_map=null)
```

Default outputs under `stats/` are EDPL, KRD, and edge-PCA. `--legacy_stats`
additionally runs the preserved ADCL, guppy EDPL/PCA/KR, and alpha-diversity tools
under `stats/legacy/`; require refpkg/map inputs only where their process signatures
need them.

Standalone `stats_entry` accepts `--dedup_jplace` and a
`--specimen_jplaces` glob/directory. It does not recreate the specimen split.

## Verify

- Default invocation creates only modern outputs.
- KRD/edge-PCA consume specimen placements; EDPL consumes deduplicated placement.
- Legacy outputs appear only with the flag and required inputs.
- Both placer outputs pass through the same stats workflow.

## Commit / PR

- Title: `feat(stats): modular modern defaults with explicit legacy metrics`
