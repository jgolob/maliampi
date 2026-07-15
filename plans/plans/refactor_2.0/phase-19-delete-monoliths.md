# Phase 19 — Delete the old monoliths

**Flag:** 🟢 SIMPLE · **Branch:** `refactor2/19-delete-monoliths` ·
**Files:** delete `modules/epang_place_classify.nf` and
`modules/pplacer_place_classify.nf` · **Depends on:** 18 and 25 merged

## Goal

Delete the monoliths only after all workflows and CI jobs use the modular entry
points. Do not coordinate deletion and CI migration in the same uncertain window;
Phase 25 must land first.

## Pre-flight

```bash
git grep -n "epang_place_classify\|pplacer_place_classify" -- '*.nf' '.github/**'
```

Proceed only when hits are confined to the files being deleted or explanatory
history/docs. Verify no process definition will disappear without a canonical or
explicit legacy replacement.

## Verify

- Repository search finds no live imports/jobs.
- Every standalone entry and the composed workflow parse.
- Phase-25 CI remains green after deletion.

## Commit / PR

- Title: `refactor: remove absorbed placement/classification monoliths`
