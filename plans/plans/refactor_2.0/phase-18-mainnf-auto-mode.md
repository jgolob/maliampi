# Phase 18 — Thin full-workflow composition

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/18-main-compose` · **Files:** `main.nf` ·
**Depends on:** 13–17

## Goal

Keep the rare FASTQ-to-terminal-output workflow, but make it a thin composition of
the same independently runnable stages.

## Flow

1. Run `sv`.
2. Use `--refpkg` after strict validation, or build and validate a new refpkg from
   explicit repository inputs.
3. Run `place` with EPA-ng by default or pplacer when selected.
4. Fan out placement into taxonomy and modern stats in parallel.
5. Phase 22 adds phylotypes as a third default-on peer branch.

Add independent `--skip_taxonomy`, `--skip_stats`, and (after Phase 22)
`--skip_phylotypes`. Selecting pplacer changes only placement. Legacy taxonomy,
reduplication, and stats remain explicit options.

When an existing refpkg is provided, gate the FASTQ manifest on successful refpkg
validation. Repository FASTA/sequence-info are required only for the build path;
email is conditional on an actual download operation.

## Verify

- Build-refpkg and provided-refpkg paths both work from the same stage interfaces.
- The provided-refpkg path does not create a new package.
- Each skip flag disables only its peer branch.
- Outputs use `sv/`, `refpkg/`, `placement/`, `taxonomy/`, and `stats/`.

## Commit / PR

- Title: `feat(main): thin modular composition around validated refpkg`
