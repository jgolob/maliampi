# Phase 05 — Tidy `pplacer_place_classify.nf` + `preprocess.nf` shebangs

**Flag:** 🟢 SIMPLE · **Branch:** `refactor2/05-tidy-pplacer-preprocess` ·
**Files:** `modules/pplacer_place_classify.nf`, `modules/preprocess.nf` ·
**Depends on:** 01

## Goal
Normalize the largest remaining batch of old-style directives (31 in pplacer) and
fix shebangs. **No behavior change.**

> Why bother with pplacer if it gets absorbed later? Because Track 3 *moves* its
> process bodies into `place/classify/stats`. Normalizing here means those moves
> are clean copies, keeping the 🔴 reorg diffs small.

## Steps
1. **`pplacer_place_classify.nf` — normalize all 31 directives:**
   ```bash
   grep -nE "^\s*(container|label|errorStrategy|publishDir|cpus|memory|tag)\s*=" modules/pplacer_place_classify.nf
   ```
   Mechanical: remove `=`, single-quote `label`/`errorStrategy`/`publishDir` string
   args. Leave the script bodies untouched.
2. **python3 shebangs in both files:**
   ```bash
   grep -n "usr/bin/env python$" modules/pplacer_place_classify.nf modules/preprocess.nf
   ```
   Change each `python` → `python3`.

## Verify
```bash
grep -nE "^\s*(container|label|errorStrategy)\s*=" modules/pplacer_place_classify.nf   # empty
grep -n "usr/bin/env python$" modules/pplacer_place_classify.nf modules/preprocess.nf  # empty
nextflow run modules/pplacer_place_classify.nf --help 2>&1 | head
```

## Done-with-Track-1 check
`modules/epang_place_classify.nf` is **intentionally excluded** — it still has
old-style directives + bare `python` shebangs, and stays that way until Track 3
moves its processes out (Phases 10-12 normalize style during the move) and deletes
it (Phase 19). Normalizing it here would be wasted/discarded work.
```bash
git grep -nE "^\s*(container|label|errorStrategy)\s*=" -- 'modules/*.nf' ':!modules/epang_place_classify.nf'   # empty
git grep -n "usr/bin/env python$" -- 'modules/*.nf' ':!modules/epang_place_classify.nf'                        # empty
```

## Commit / PR
- Title: `style(pplacer,preprocess): normalize directives + python3 shebangs`
