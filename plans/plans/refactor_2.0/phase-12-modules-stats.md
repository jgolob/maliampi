# Phase 12 — Modern and legacy stats processes

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/12-modules-stats` ·
**Files:** `modules/stats.nf`, both monolith includes · **Depends on:** 11

## Goal

Centralize placement statistics while distinguishing modern default gappa metrics
from preserved guppy/pplacer compatibility metrics.

## Processes

Modern default group:

- `EDPL`
- `Gappa_KRD`
- `Gappa_ePCA`

Legacy opt-in group:

- `PplacerADCL`
- `PplacerEDPL`
- `PplacerPCA`
- `PplacerAlphaDiversity`
- `CombineSpAd`
- `PplacerKR`

## Steps

1. Move and deduplicate process definitions without changing commands.
2. Mark the two groups in the module but defer default/opt-in orchestration and the
   clean `stats/` publication layout to Phase 15.
3. Preserve exact process inputs: deduplicated metrics use `dedup.jplace`; KRD,
   edge-PCA, and per-specimen alpha calculations use specimen jplace files.

## Verify

- One definition per process.
- Modern and legacy fixture outputs retain their schemas and expected numeric
  tolerances.
- No process silently substitutes deduplicated placement for specimen placement.

## Commit / PR

- Title: `refactor(stats): centralize modern and legacy placement metrics`
