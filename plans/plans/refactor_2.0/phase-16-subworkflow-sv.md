# Phase 16 — Standalone FASTQ-to-SV subworkflow

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/16-subwf-sv` ·
**Files:** `subworkflows/local/sv.nf` · **Depends on:** 06, 12a

## Goal

Turn a manifest into the canonical four-file SV contract while preserving manifest,
preprocessing, DADA2, and failure-reporting behavior.

## Contract

`sv(manifest)` and `sv_entry --manifest` emit:

- `sv/sv.fasta`
- `sv/sv.map.csv`
- `sv/sv.weights.csv`
- `sv/sv.long.csv`
- The filtered seqtable and failure report as secondary outputs

Use the Phase-12a normalizer rather than publishing DADA2-specific names as the
public stage API. Accept optional `--dataset_id` to namespace SV identifiers when
the output will enter a cumulative phylotype set.

## Steps

1. Move `read_manifest`, preprocessing, DADA2, and `output_failed` wiring from
   `main.nf` without changing failure semantics.
2. Normalize DADA2 outputs to the canonical files and validate cross-file IDs and
   counts.
3. Keep `read_manifest` as a function; do not pretend it is a Nextflow process.

## Verify

- Fixture FASTQs produce the canonical four files.
- Counts and specimen failure classifications match existing fixtures.
- A dataset namespace relabels all four files consistently.

## Commit / PR

- Title: `feat(sv): standalone canonical FASTQ-to-SV stage`
