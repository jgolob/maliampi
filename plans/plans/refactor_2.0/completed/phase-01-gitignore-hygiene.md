# Phase 01 — gitignore + remove stray runtime artifacts

**Flag:** 🟢 SIMPLE · **Branch:** `refactor2/01-gitignore` · **Files:** `.gitignore`
· **Depends on:** none

## Goal
Stop Nextflow runtime junk from polluting the repo and clear what's already there.

## Background
`.gitignore` currently only has `*working/` and `*Audit` for Nextflow-ish paths.
The working tree has untracked `.nextflow/`, `.nextflow.log.1`..`.5`, `work/`,
`output/`.

## Steps
1. Append a Nextflow block to `.gitignore`:
   ```gitignore
   # Nextflow runtime
   .nextflow/
   .nextflow.log*
   .nextflow.pid
   work/
   output/
   results/
   ```
2. Remove the stray untracked artifacts (they are regenerable run output — confirm
   none are intentional fixtures, which they are not; fixtures live under `tests/`):
   ```bash
   rm -rf .nextflow .nextflow.log.* work output
   ```
   > Do **not** delete `tests/` or `maliampi-practice-data/`.

## Verify
```bash
git status --porcelain        # no .nextflow*/work/output noise
git check-ignore work output .nextflow   # all three printed = ignored
```

## Commit / PR
- Title: `chore: gitignore nextflow runtime artifacts`
- Body: one line; mention removed stray `work/`, `output/`, `.nextflow*`.
