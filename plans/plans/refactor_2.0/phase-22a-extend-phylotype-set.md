# Phase 22a — Validate and extend a phylotype set from new FASTQs

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/22a-extend-set` ·
**Files:** `modules/phylotypes.nf`, `subworkflows/local/extend_phylotype_set.nf` ·
**Depends on:** 09a, 13–16, 21, 22

## Goal

Provide the routine workflow for new FASTQs against an existing refpkg and
phylotype set. Validate all compatibility before FASTQ work, place the new SVs,
extend every selected mapping without unassigned SVs, and emit a new versioned set.

## Standalone contract

`extend_phylotype_set_entry` requires:

- `--manifest`
- `--dataset_id` (namespaces new SV IDs and prevents cumulative collisions)
- `--refpkg`
- `--phylotype_set`
- Optional `--placer epang|pplacer` and threshold additions/overrides

It may also run taxonomy and stats for the new data, controlled by the same peer
skip flags as main.

## Validate before work

`ValidatePhylotypeSet` must verify:

- Set-manifest checksums and content-derived `set_id`.
- Set `refpkg_id` equals the validated package ID.
- Cumulative jplace tree fingerprint matches the refpkg/placement contract.
- Each stored threshold mapping covers cumulative jplace SVs exactly once.
- Mapping parameters and checksums match the manifest.
- New dataset/SV identifiers do not collide with cumulative identifiers.

Gate the FASTQ manifest on the combined successful refpkg+set validation token.

## Extension behavior

1. Run canonical SV creation and placement against the validated refpkg.
2. Merge old and new jplace records only after tree/field normalization and duplicate
   ID checks; retain one cumulative jplace.
3. For every threshold already in the set, run `extend_phylotypes` and preserve all
   existing IDs while assigning every new SV to an old or newly created group.
4. For a newly requested threshold, run full phylotype generation over the merged
   cumulative jplace to initialize that threshold.
5. Keep unrequested stored thresholds; never delete them implicitly.
6. Emit new-run mapping/count/relative-abundance tables for every selected threshold.
7. Emit a new immutable set bundle with `parent_set_id`, new `set_id`, cumulative
   jplace, complete mappings, checksums, and provenance. Never overwrite the parent.

## Invariants

- Every cumulative SV appears exactly once at every stored threshold.
- Every new read count is represented in the new-run tables.
- There is no orphan/unassigned state and no silent row loss.
- Existing phylotype IDs and memberships are stable except for adding new members.

## Verify

Test join-existing, create-new, add-new-threshold, repeat-extension, refpkg mismatch,
tree mismatch, identifier collision, tampered bundle, and count conservation.

## Commit / PR

- Title: `feat(phylotypes): validated versioned-set extension from new FASTQs`
