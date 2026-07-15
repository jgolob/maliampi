# Phase 10 — Placement processes and engine adapters

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/10-modules-place` ·
**Files:** `modules/place.nf`, both placement/classification monoliths ·
**Depends on:** 09

## Goal

Move shared alignment, EPA-ng placement, pplacer placement, and specimen splitting
into one placement module. EPA-ng and pplacer are interchangeable placement engines,
not selectors for different taxonomy or stats stacks.

## Processes

Move one canonical definition of:

- `AlignSV`
- `CombineAln_SV_refpkg`
- `ConvertAlnToFasta`
- `EPAngPlacement`
- `PplacerPlacement`
- `MakeSplit`
- `GappaSplit`

Leave SV table adapters and `PplacerReduplicate` for Phase 12a's conversion module.
Leave taxonomy and stats for Phases 11–12.

## Steps

1. Copy process bodies from both monoliths, reconciling duplicates against their
   actual signatures. Do not preserve two copies merely because both files had one.
2. Normalize Python shebangs only; scientific command changes are out of scope.
3. Keep the old workflows temporarily wired to the new definitions.
4. Do not add backend-specific downstream calls to `modules/place.nf`.

## Verify

- Every moved process has exactly one definition.
- Existing EPA-ng and pplacer fixture commands still emit valid deduplicated
  jplace files with the original SV identifiers.
- The same refpkg produces a compatible jplace tree fingerprint from both engines;
  placement likelihoods are not expected to be numerically identical.

## Commit / PR

- Title: `refactor(place): centralize EPA-ng and pplacer placement processes`
