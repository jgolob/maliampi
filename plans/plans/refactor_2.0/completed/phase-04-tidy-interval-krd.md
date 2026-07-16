# Phase 04 — Tidy `modules/interval_krd.nf` + legacy banner

**Flag:** 🟢 SIMPLE · **Branch:** `refactor2/04-tidy-interval-krd` ·
**Files:** `modules/interval_krd.nf` · **Depends on:** 01

## Goal
Normalize the 4 old-style directives and mark legacy. **No behavior change.**

## Steps
1. **Normalize directives:**
   ```bash
   grep -nE "^\s*(container|label|errorStrategy|publishDir|cpus|memory|tag)\s*=" modules/interval_krd.nf
   ```
   Delete `=`, single-quote `label`/`errorStrategy`.
2. **Legacy banner** under the `dsl=2` line:
   ```nextflow
   // LEGACY — standalone interval KR-distance utility. Not wired into the
   // 5-stage flow. Kept for reference; unmaintained.
   ```
   (This file has no python shebangs to change.)

## Verify
```bash
grep -nE "^\s*(container|label|errorStrategy)\s*=" modules/interval_krd.nf   # empty
nextflow run modules/interval_krd.nf --help 2>&1 | head    # parses
```

## Commit / PR
- Title: `style(interval_krd): normalize directives, mark legacy`
