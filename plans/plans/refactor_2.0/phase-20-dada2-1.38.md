# Phase 20 — DADA2 → 1.38 + validate Rscripts + refresh fixtures

**Flag:** 🔴 MAJOR (likely debugging + fixture regen) · **Branch:**
`refactor2/20-dada2-138` · **Files:** `modules/dada2.nf`, `tests/sv/*` ·
**Depends on:** Track 1 merged (independent of the reorg; can run in parallel)

## Goal
Modernize DADA2 from 1.26.0 → the current latest stable biocontainer and prove the
SV outputs are still correct.

## Target image
```
params.container__dada2 = "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"
```
(Verify it's still the latest at implementation:
`curl -s 'https://quay.io/api/v1/repository/biocontainers/bioconductor-dada2/tag/?limit=10&onlyActiveTags=true'`.)

## Steps
1. Bump `params.container__dada2` (and any related dada2 container default) in
   `modules/dada2.nf`.
2. **Validate the Rscripts.** These calls are API-stable across 1.26→1.38 but
   re-check argument names against the 1.38 docs: `filterAndTrim`, `derepFastq`,
   `learnErrors`, `dada`, `mergePairs`, `makeSequenceTable`, `removeBimeraDenovo`.
   Watch for: changed defaults, deprecation warnings treated as errors, R 4.5 vs 4.2
   behavior.
3. **Bridge-image risk:** `golob/dada2-pplacer` and `golob/dada2-fast-combineseqtab`
   consume DADA2's seqtab/RDS. Confirm they still parse 1.38 output. If the
   serialized format shifted, rebuild those images (their Dockerfiles should be in
   `docker/`; if not, this phase grows — flag it).
4. **Refresh fixtures.** A minor-version DADA2 bump can legitimately change SVs.
   Run the dada2 workflow on practice data, inspect the diff vs `tests/sv/`, and if
   the changes are explainable (not a bug), update the fixtures:
   ```bash
   nextflow run subworkflows/local/sv.nf -entry sv_entry -profile test,docker \
     --manifest maliampi-practice-data/manifest.csv --output /tmp/sv_new
   # eyeball /tmp/sv_new/sv vs tests/sv ; if good, copy over and commit
   ```

## Verify
- dada2 workflow runs clean on practice data, no R errors/warnings-as-errors.
- Downstream placement and taxonomy still consume the new SVs (run the relevant
  modular entries; full composition is secondary).
- Fixture diff reviewed and justified in the PR body.

## Commit / PR
- Title: `feat(dada2): upgrade to bioconductor-dada2 1.38 + refresh fixtures`
- 🔴 Body MUST explain any fixture changes (why the SVs differ).
