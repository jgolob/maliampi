# Refactor 2.0 — Execution Plan (granular, phase-by-phase)

This directory is currently **workspace-local and untracked**. Decide explicitly
whether to commit or exclude it before any Git operation. It holds one task file
per phase. Each phase ends with its own reviewed PR. Do not commit, push, publish an
image, or start a later phase without explicit approval.

## Target architecture

Modular entry points are the primary product. `main.nf` is only a thin, rarely
used composition.

```text
FASTQ -> canonical SV files -> validated refpkg -> placement
                                               |-> taxonomy
                                               |-> phylotypes (many thresholds)
                                               `-> stats
```

- The refpkg is the immutable compatibility anchor and has a content-derived
  `refpkg_id`.
- EPA-ng is the default placement engine; pplacer is a fallback placement engine.
  Both emit the same deduplicated and per-specimen jplace contract.
- Taxonomy, phylotypes, and stats are independent peer pathways.
- Phylotype sets are validated, versioned, and extensible under one refpkg. New
  SVs either join an existing phylotype or form a new one; none are orphaned,
  unassigned, or dropped.
- Nextflow is pinned to **26.04.6** in local instructions, CI, and nf-test.

## How to use this

1. Read this README, then the next phase file and every dependency it names.
2. Create only the phase branch and change only the files in scope.
3. Run the phase's verification after receiving approval for any resource-heavy or
   externally mutating action.
4. Open a PR and wait for review/merge unless explicitly told to stack work.
5. Record extra cleanup for a later phase instead of expanding the active PR.

## Conventions

- Process directives use space form and single-quoted labels.
- Container parameters are `params.container__tool`, declared once in committed
  configuration rather than copied into every module.
- Scripts use `#!/usr/bin/env python3`.
- Container images are pinned to immutable versions/builds; never use `latest`.
- Imported workflows take explicit channels and emit named channels. Only
  standalone entry wrappers translate `params` into channels.
- Branches use `refactor2/NN-short-slug` (including `09a`, `12a`, etc.).
- A behavior-preserving move must not silently become a scientific change.

## Phase index

### Track 1 — Completed hygiene

| Phase | Title | Flag |
|---|---|---|
| 01–05a | Git hygiene, module cleanup, Nextflow 26 compatibility | complete |

### Track 2 — Config and lint scaffolding

| Phase | Title | Flag |
|---|---|---|
| 06 | Committed config, centralized defaults, Nextflow 26.04.6 | 🟡 |
| 07 | EditorConfig + Prettier | 🟢 |
| 08 | nf-core lint scaffolding | 🟡 |

### Track 3 — Modular reorganization

| Phase | Title | Flag |
|---|---|---|
| 09 | Shared refpkg extraction utility | 🟡 |
| 09a | Refpkg identity + validation entry | 🔴 |
| 10 | Placement processes and engine adapters | 🔴 |
| 11 | Taxonomy module + isolated legacy pplacer taxonomy | 🔴 |
| 12 | Modern and legacy stats processes | 🔴 |
| 12a | Canonical conversion module + named entries | 🔴 |
| 13 | Standalone placement subworkflow | 🔴 |
| 14 | Standalone taxonomy subworkflow | 🟡 |
| 15 | Standalone modern/legacy stats subworkflow | 🟡 |
| 16 | Standalone FASTQ-to-SV subworkflow | 🟡 |
| 17 | Standalone build-and-validate-refpkg subworkflow | 🟡 |
| 18 | Thin full-workflow composition | 🔴 |
| 19 | Delete old monoliths, only after Phase 25 | 🟢 |

### Track 4 — Feature work

| Phase | Title | Flag |
|---|---|---|
| 20 | DADA2 1.38 + fixture refresh | 🔴 |
| 20a | Upstream phylotypes extension API + v2.1.0 release | 🔴 |
| 21 | Pinned phylotypes container | 🟡 |
| 22 | Parallel multi-threshold phylotypes pathway | 🔴 |
| 22a | Validate and extend a versioned phylotype set from new FASTQs | 🔴 |
| 23 | Biocontainer sweep + interval-KRD gappa migration | 🟡/🔴 |

### Track 5 — Tests

| Phase | Title | Flag |
|---|---|---|
| 24 | Content, identity, and conservation assertions | 🟡 |
| 25 | Core per-entry CI; migrate old monolith jobs | 🔴 |
| 25a | Phylotypes and set-extension entry CI | 🔴 |
| 26 | Seed nf-test with functions, stages, and set extension | 🔴 |

### Track 6 — Docs and schema

| Phase | Title | Flag |
|---|---|---|
| 27 | Container inventory and provenance | 🟢 |
| 28 | Changelog, contributing, modular/refpkg/set docs | 🟢 |
| 29 | Conditional parameter schema | 🟡 |

## Required ordering

Run Tracks 1→2→3. Add Phases 09a and 12a before Phase 13. Run Phase 25 after
Phase 18 and **before Phase 19**, then delete the monoliths. Phase 20a must release
the upstream extension API before Phase 21 pins it; Phase 22a depends on 09a, 16,
21, and 22, and Phase 25a adds its CI coverage. Documentation and schema record the
final interfaces last.
