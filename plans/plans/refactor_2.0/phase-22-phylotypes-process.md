# Phase 22 — Parallel multi-threshold phylotypes pathway

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/22-phylotypes-pathway` ·
**Files:** `modules/phylotypes.nf`, `subworkflows/local/phylotypes.nf`, `main.nf` ·
**Depends on:** 13, 16, 21

## Goal

Create phylotypes as a first-class peer to taxonomy. For every configured threshold,
emit the SV-to-phylotype mapping and taxonomy-like specimen count/relative-abundance
tables. Also package the initial versioned phylotype set.

## Parameters

```nextflow
params.phylotype_thresholds = [0.1, 0.3, 1.0] // configurable defaults
params.phylotype_add_thresholds = []           // append without restating defaults
params.phylotype_lwr_overlap = 0.1
params.phylotype_distance = 'legacy'
```

Normalize the effective threshold list by parsing numbers, requiring positive
finite values, deduplicating, and sorting. An explicit `phylotype_thresholds`
replaces the defaults; `phylotype_add_thresholds` appends to the selected/base list.
KR is opt-in and requires an explicitly supplied calibrated threshold list.

## Processes and contract

`phylotypes(dedup_jplace, sv_long, placement_identity, dataset_id, thresholds, ...)`
scatters tuples `(threshold, dedup_jplace)` into:

1. `Phylotypes`: complete `phylotype,sv` mapping for one threshold.
2. `MakeWidePhylotypeTables`: join the mapping to `sv_long`; emit per-specimen
   read counts and relative abundances.
3. `BuildPhylotypeSet`: package cumulative jplace, mappings, and a manifest bound
   to `refpkg_id`.

Publish per threshold:

```text
phylotypes/threshold-0p1/sv_phylotype.csv
phylotypes/threshold-0p1/phylotype_wide_nreads.csv
phylotypes/threshold-0p1/phylotype_wide_ra.csv
```

Repeat for every effective threshold. The initial set manifest records `set_id`,
`parent_set_id=null`, `refpkg_id`, dataset ID, cumulative-jplace checksum/tree
fingerprint, phylotypes version and parameters, thresholds, and mapping checksums.
Compute `set_id` as SHA-256 of the canonical manifest content excluding `set_id`
itself, with all referenced artifact SHA-256 values included and keys serialized in
stable order.

Add phylotypes to `main.nf` as a default-on peer branch controlled only by
`--skip_phylotypes`. It must not live in or depend on the taxonomy module.

## Verify

- Every jplace SV occurs exactly once in every threshold mapping.
- Count tables conserve each specimen's `sv_long` total.
- Every nonempty relative-abundance row sums to 1 within tolerance.
- Output discovery is driven by the effective threshold list, not three hardcoded
  paths.
- Main runs taxonomy, phylotypes, and stats in parallel by default.

## Commit / PR

- Title: `feat(phylotypes): parallel multi-threshold mappings and abundance tables`
