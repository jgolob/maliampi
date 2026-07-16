# Phase 11 — Taxonomy module + isolated legacy taxonomy

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/11-modules-taxonomy` ·
**Files:** `modules/taxonomy.nf`, `modules/legacy/pplacer_taxonomy.nf`, old
monolith includes · **Depends on:** 10

## Goal

Make gappa taxonomy the common supported taxonomy pathway for placements from
either engine. Preserve the old pplacer classifier only as an explicit EoL backup;
`--placer pplacer` must never select it automatically.

## Processes

Move to `modules/taxonomy.nf`:

- `MakeEPAngTaxonomy`
- `Gappa_Classify`
- `Gappa_Extract_Taxonomy`
- `Make_Wide_Tax_Table`
- `Extract_Taxonomy`

Move to `modules/legacy/pplacer_taxonomy.nf`:

- `ClassifyDB_Prep`
- `ClassifySV`
- `ClassifyMCC`
- `ClassifyTables`

## Steps

1. Move process bodies from both monoliths and remove duplicate definitions.
2. Keep command behavior and existing publication paths during this pure move;
   Phase 14 intentionally changes the public directory to `taxonomy/`.
3. Add a clearly marked standalone legacy wrapper only if needed to preserve a
   currently exercised entry point. Do not include legacy taxonomy from main.
4. Do not add phylotypes here; it is a peer pathway in Phase 22.

## Verify

- One definition per process.
- Gappa taxonomy accepts fixture jplace from both placers.
- Existing taxonomy fixture content matches after stable sorting.
- Choosing pplacer alone does not invoke any legacy taxonomy process.

## Commit / PR

- Title: `refactor(taxonomy): common gappa path and isolated legacy classifier`
