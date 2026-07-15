# Phase 05a (unplanned) — Nextflow 26.x compatibility sweep

**Flag:** retroactive record · **Branch:** `refactor2/01-track1-hygiene` (same PR as
01-05) · **Files:** `main.nf`, all of `modules/*.nf`, `.github/workflows/test.yaml`
· **Depends on:** 01-05 · **Status:** ✅ DONE (merged in PR #19, commits
`5a7fcf9..ebcbec9`)

## Why this exists
Not in the original plan. CI installs Nextflow via
`wget -qO- get.nextflow.io | bash` (unpinned), which started pulling **26.04.3**
mid-refactor. That version enforces much stricter DSL2 parsing than this codebase
was written against, and the Track 1 PR stopped compiling in CI. Fixing it touched
every module file, so it's recorded here rather than silently folded into Phase 05.

## What changed (all behavior-preserving)
- **CI runner/toolchain**: bumped `ubuntu-20.04` → `ubuntu-22.04`, added explicit
  Java 17 setup before the Nextflow install, made the docker-image cleanup step
  tolerate an empty image list.
- **`include { x } from '...' params (...)`** (deprecated/removed syntax) — removed
  entirely from `main.nf`. Each module already self-declares its own `params.*`
  defaults and reads the global `params` scope directly, so the include blocks were
  redundant.
- **Workflow section syntax**: repeated `input: name` / `take: name` lines →
  single `take:` label followed by bare identifiers, one per line (required by
  26.x). Fixed in `dada2.nf`, `preprocess.nf`, `swarm.nf`, `refpackage.nf`,
  `pplacer_place_classify.nf`, `epang_place_classify.nf`.
- **`container__xxx = "image:tag"` → `params.container__xxx = "image:tag"`**
  (and all ~100 `${container__xxx}` usages → `${params.container__xxx}`). Top-level
  bare variable assignments are no longer legal when mixed with process/workflow
  declarations; `params.x = value` is the only allowed form. **This is now the
  convention for every container pin** — see updated README Conventions.
- **`script:` labels**: every process body (`"""..."""`) preceded by `output:` (or
  other sections) now has an explicit `script:` label, as 26.x requires.
- **Old-style directives** (`container =`, `label =`, `errorStrategy =`, `cache =`,
  `afterScript "..."`) normalized to space-form, single-quoted — including in
  `epang_place_classify.nf`, which Phase 05 had explicitly deferred to Track 3.
  **`epang_place_classify.nf`'s directives are therefore already clean**; only its
  10 bare `#!/usr/bin/env python` shebangs remain for Phases 10-12 to fix during the
  move.
- **`dada2.nf` closure/channel fixes**: `def` added to closure-local vars; the old
  `.tap { a; b }` implicit-variable pattern (no longer supported) replaced with
  `.set { dada2_demultiplex_dada_split }` reused twice (DSL2 channels are
  re-consumable); two `publishDir "...${batch}"` directives wrapped in closures
  (`publishDir { "...${batch}" }`) since input-var interpolation in directives is
  now lazy-only.
- **`afterScript 'rm -r refpkg/'` / `'rm -rf dl/'`** (10 occurrences across
  `refpackage.nf`, `epang_place_classify.nf`, `pplacer_place_classify.nf`) → appended
  `|| true`. Docker tasks write `refpkg/`/`dl/` as root; the host-side afterScript
  cleanup as the unprivileged runner failed with "Permission denied" and was
  failing otherwise-successful tasks.

## Verification performed
`nextflow run <file> -preview` clean (no `Error <file>:<line>` blocks) for every
module plus `main.nf`, in dependency order: `dada2.nf` → `swarm.nf` →
`refpackage.nf` → `preprocess.nf` → `interval_krd.nf` → `pplacer_place_classify.nf`
→ `epang_place_classify.nf` → `main.nf`. Full `pplacer_place_classify.nf` module run
end-to-end with `-with-docker` against `tests/refpkg/refpkg.tar.gz` +
`tests/sv/*` fixtures: 20/20 processes succeeded.
(`interval_krd.nf`'s one remaining `-preview` error,
`Missing process or function plus([/*.jplace.gz])`, is pre-existing —
`params.old_jplace` defaults to `false`; expected unless real
`--old_jplace`/`--new_jplace` are passed.)

## Follow-up (not done here, candidates for Phase 06)
- CI still installs Nextflow unpinned (`wget -qO- get.nextflow.io | bash`). Pinning
  `NXF_VER` would prevent a future Nextflow release from breaking the pipeline again
  with no warning — worth adding alongside the Phase 06 config work.
