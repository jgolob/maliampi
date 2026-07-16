# Phase 29 — Conditional Nextflow parameter schema

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/29-schema` ·
**Files:** `nextflow_schema.json` · **Depends on:** 22a, 28

## Goal

Encode the modular v2 interfaces and their conditional requirements rather than
maintaining a flat list of nominally required params.

## Groups

- Input/output and dataset namespace
- Canonical SV conversion
- Reference package build/validation
- Placement
- Taxonomy
- Phylotypes and versioned-set extension
- Stats and explicit legacy options

## Required validation rules

- Use `oneOf` for canonical SV source modes.
- Require either an existing refpkg or the complete repository inputs needed to
  build one; do not require both.
- Require email/credentials only for the download path that uses them.
- Restrict placer to `epang|pplacer`; it does not select taxonomy/stats behavior.
- Define `phylotype_thresholds` and `phylotype_add_thresholds` as arrays of positive
  finite numbers. Defaults are `[0.1, 0.3, 1.0]`, but remain replaceable/appendable.
- Require an explicit calibrated threshold list when distance is `kr`.
- Require refpkg and map only for legacy stats that consume them.
- Require `manifest`, `dataset_id`, `refpkg`, and `phylotype_set` for set extension.
- Reject contradictory skip/run/legacy combinations and unsafe dataset IDs.

## Verify

```bash
nf-core pipelines schema lint
```

Add positive and negative examples for every `oneOf`/conditional branch, especially
threshold replacement/addition and phylotype-set extension.

## Commit / PR

- Title: `feat(schema): modular conditional and set-extension parameters`
