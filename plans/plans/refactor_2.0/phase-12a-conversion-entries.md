# Phase 12a — Canonical conversion module and named entries

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/12a-conversion` ·
**Files:** `modules/conversion.nf`, `subworkflows/local/convert.nf` ·
**Depends on:** 09a, 10

## Goal

Create one compatibility boundary for legacy SV tables, placement alignment
formats, and pplacer reduplication. Downstream stages accept only canonical files.

## Named entries

`normalize_sv_entry` accepts exactly one source form:

1. `--seqtable`
2. `--sv_fasta --sharetable`
3. `--sv_fasta --sv_map --sv_weights`
4. `--sv_fasta --sv_long`
5. All four canonical files as a validated pass-through

It emits `sv/sv.fasta`, `sv/sv.map.csv`, `sv/sv.weights.csv`, and
`sv/sv.long.csv`. Require identifier/count consistency and reject mixed modes.

`prepare_alignment_entry --sv_fasta --refpkg` validates the package, then emits
`alignment/combined.sto` and `alignment/combined.fasta` through shared refpkg
extraction and alignment processes.

`pplacer_reduplicate_entry --dedup_jplace --sv_weights` emits the explicit legacy
artifact `placement/redup.jplace.gz`.

## Dataset namespace

Accept `--dataset_id` when creating data intended for a cumulative phylotype set.
Relabel FASTA, map, weights, and long-table SV identifiers consistently so
independent runs cannot collide. Reject an unsafe/empty namespace.

## Verify

- Every accepted source form yields the same canonical schema and identifiers.
- Invalid/mixed forms fail before a converter starts.
- FASTA IDs equal the long-table `sv` values and representative weight IDs.
- Alignment FASTA and Stockholm outputs contain the same records.
- Reduplication remains independently runnable and is not required by EPA-ng.

## Commit / PR

- Title: `feat(conversion): canonical SV, alignment, and legacy redup entries`
