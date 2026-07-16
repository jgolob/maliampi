# Phase 03 — Tidy `modules/swarm.nf` + legacy banner

**Flag:** 🟢 SIMPLE · **Branch:** `refactor2/03-tidy-swarm` ·
**Files:** `modules/swarm.nf` · **Depends on:** 01

## Goal
Normalize directives + shebangs, and clearly mark this module **legacy** (swarm is
an alternative ASV method, not part of the supported 5-stage flow). **No behavior
change.** (The `golob/gappa` → biocontainer migration is a *separate* phase, 23 —
do **not** change the container here.)

## Steps
1. **Normalize the 8 old-style directives** (delete `=`, single-quote
   `label`/`errorStrategy`):
   ```bash
   grep -nE "^\s*(container|label|errorStrategy|publishDir|cpus|memory|tag)\s*=" modules/swarm.nf
   ```
2. **python3 shebangs:**
   ```bash
   grep -n "usr/bin/env python$" modules/swarm.nf   # change each to python3
   ```
3. **Add a legacy banner** as the first comment block under the shebang/`dsl=2`
   line:
   ```nextflow
   // LEGACY — alternative ASV method (swarm). Not part of the supported
   // 5-stage flow (sv→place→classify→stats). Kept for reference; unmaintained.
   ```

## Verify
```bash
grep -nE "^\s*(container|label|errorStrategy)\s*=" modules/swarm.nf   # empty
grep -n "usr/bin/env python$" modules/swarm.nf                       # empty
nextflow run modules/swarm.nf --help
```

## Commit / PR
- Title: `style(swarm): normalize DSL2 directives, python3, mark legacy`
