# Phase 02 — Tidy `modules/refpackage.nf`

**Flag:** 🟢 SIMPLE · **Branch:** `refactor2/02-tidy-refpackage` ·
**Files:** `modules/refpackage.nf` · **Depends on:** 01

## Goal
Finish the DSL2 directive normalization in this file and drop migration cruft.
**No behavior change.**

## Steps
1. **Normalize the 2 remaining old-style directives** to space-form, single quotes:
   ```bash
   grep -nE "^\s*(container|label|errorStrategy|publishDir|cpus|memory|tag)\s*=" modules/refpackage.nf
   ```
   For each hit, delete the `=` (e.g. `errorStrategy = 'finish'` →
   `errorStrategy 'finish'`; `container = "${x}"` → `container "${x}"`). Convert any
   double quotes on `errorStrategy`/`label` to single quotes.
2. **Remove the commented old container line** near the top (the
   `//container__raxmlng = '…1.0.3…'` line that sits above the active 1.2.2 pin).
   Keep the active line.
3. **python3 shebangs:**
   ```bash
   grep -n "usr/bin/env python$" modules/refpackage.nf
   ```
   Change each `#!/usr/bin/env python` → `#!/usr/bin/env python3`.

## Verify
```bash
grep -nE "^\s*(container|label|errorStrategy)\s*=" modules/refpackage.nf   # empty
grep -n "usr/bin/env python$" modules/refpackage.nf                       # empty
nextflow run modules/refpackage.nf --help    # parses (prints help), no DSL errors
```

## Commit / PR
- Title: `style(refpackage): normalize DSL2 directives + python3 shebangs`
