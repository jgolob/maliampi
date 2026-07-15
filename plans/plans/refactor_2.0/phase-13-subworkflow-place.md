# Phase 13 — Standalone placement subworkflow

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/13-subwf-place` ·
**Files:** `subworkflows/local/place.nf`, `modules/place.nf` as required ·
**Depends on:** 09a, 10, 12a

## Goal

Expose EPA-ng and pplacer as interchangeable placement engines behind one modular
contract. Everything after the placement command is common.

## Contract

Imported workflow inputs are explicit:

```text
place(sv_fasta, sv_long, validated_refpkg, placer,
      legacy_redup=false, sv_weights=null)
```

Standalone `place_entry` accepts canonical `--sv_fasta`, `--sv_long`, `--refpkg`,
and optional `--placer epang|pplacer`. It validates the refpkg internally. Legacy
table forms belong to Phase 12a and are not accepted here.

Emits:

- `dedup_jplace`
- `specimen_jplaces`
- `placement_identity` containing `refpkg_id`, tree fingerprint, placer, and tool
  version
- Optional `redup_jplace` only when `--legacy_redup` and weights are provided

## Wiring

1. Prepare the shared alignment once, producing combined Stockholm and FASTA.
2. Dispatch only the placement command:
   - EPA-ng consumes the validated refpkg components and combined FASTA.
   - pplacer consumes the combined Stockholm and refpkg archive.
3. Run common `MakeSplit(sv_long)` and `GappaSplit(dedup_jplace, split_file)`.
4. Write all engine-independent placement outputs under `placement/`.
5. Do not invoke taxonomy, phylotypes, stats, or pplacer-specific classification.

## Verify

- Both engines emit valid deduplicated jplace with all input SVs.
- Both emit the expected specimen file set.
- Placement identity matches the validated refpkg.
- `--legacy_redup` is absent by default and works independently when requested.
- Do not require numerical equality between placement engines.

## Commit / PR

- Title: `feat(place): interchangeable placement engines and common outputs`
