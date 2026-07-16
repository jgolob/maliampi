# Phase 09 — Shared refpkg extraction utility

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/09-refpkg-util` ·
**Files:** `modules/refpkg_utils.nf`, both placement/classification monoliths ·
**Depends on:** 06

## Goal

Create one `ExtractRefpkg` definition used by every placer, taxonomy workflow, and
validator. Preserve supported refpkg formats while fixing the broken single-format
alignment branches.

## Steps

1. Move one canonical `ExtractRefpkg` process to `modules/refpkg_utils.nf` and
   remove both monolith definitions.
2. Preserve named outputs: tree, FASTA alignment, Stockholm alignment, model,
   sequence information, taxonomy, and covariance model.
3. Deliberately fix the existing defects while moving:
   - Import `Bio.AlignIO`.
   - Replace the `sopen(...)` typo with `open(...)`.
   - Support refpkgs containing both alignment formats or either one alone.
   - Reject a refpkg containing neither alignment format.
4. Use safe member lookup: reject duplicate basenames, path traversal, missing
   `CONTENTS.json`, and missing logical components.
5. Include the shared process from both old monoliths until Phase 19 deletes them.
   Do not copy shared defaults into the new module.

## Verify

```bash
rg -n "process ExtractRefpkg" modules
# exactly one definition
```

Add fixtures for both-alignment, FASTA-only, Stockholm-only, and malformed
packages. Extracted tree/alignment/model/taxonomy files must match their source
components.

## Commit / PR

- Title: `refactor(refpkg): centralize and repair extraction`
