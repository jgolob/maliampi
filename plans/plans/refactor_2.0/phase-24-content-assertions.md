# Phase 24 — Content assertions in existing CI jobs

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/24-content-assertions` ·
**Files:** `.github/workflows/test.yaml` (+ maybe `tests/assert_outputs.sh`) ·
**Depends on:** 06 and the outputs being asserted; add set assertions after 22a

## Goal
Upgrade CI from "file exists" (`[[ -s f ]]`) to "file is correct" by comparing
produced outputs against the `tests/` fixtures that already exist but are unused.

## Strategy by output type
- **Deterministic text** (`sv.map.csv`, `sv_taxonomy.csv`): exact compare after
  sorting to neutralize row order:
  ```bash
  diff <(sort output/sv/dada2.sv.map.csv) <(sort tests/sv/dada2.sv.map.csv)
  ```
- **Tables with floats / nondeterministic order** (wide tax tables, distance
  matrices): compare header + row count + a few spot values, not byte-exact.
- **jplace**: parse with `jq`/python — assert `tree` present, `placements` length
  matches, fields present. Do not byte-compare (LWR floats wobble).
- **Identity manifests**: recompute `refpkg_id`/`set_id`, component checksums, and
  tree fingerprints; do not merely assert the JSON exists.
- **Taxonomy/phylotype abundance tables**: assert specimen count conservation and
  relative-abundance row sums. Discover phylotype outputs from the configured
  threshold list rather than hardcoded filenames.

## Steps
1. Add a small `tests/assert_outputs.sh` with helpers (`assert_sorted_equal`,
   `assert_rowcount`, `assert_jplace_ok`) so each CI job calls one line.
2. In each job in `test.yaml`, replace the block of `[[ -s … ]]` lines with the
   matching assert helpers. Keep existence checks for genuinely nondeterministic
   blobs.
3. Add helpers for complete mapping partitions: every jplace SV appears once and
   only once in each threshold mapping, including after set extension.

## Verify
- CI green on a clean run.
- **Mutation test:** locally corrupt one output (e.g. `echo x >> sv.map.csv`) and
  confirm the assertion fails — proves the check bites.

## Commit / PR
- Title: `test(ci): assert output content against tests/ fixtures`
