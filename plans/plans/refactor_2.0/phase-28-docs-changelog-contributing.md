# Phase 28 — Architecture, modular, refpkg, and set documentation

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/28-docs` ·
**Files:** `CHANGELOG.md`, `CONTRIBUTING.md`, `README.md`, focused docs as needed ·
**Depends on:** 22a, 25

## Goal

Document the modular-first architecture and the two routine workflows: creating an
initial analysis under a refpkg and extending a validated phylotype set with new
FASTQs.

## Content

1. Changelog: clean v2 output API, validated/content-addressed refpkg, interchangeable
   placers, parallel taxonomy/phylotypes/stats, and versioned set extension.
2. Contributing: phase-sized PRs, directive/shebang rules, Nextflow 26.04.6,
   intentional fixture updates, and explicit approval for pushes/publication.
3. README modular table covering validation, conversion, SV, refpkg, placement,
   taxonomy, phylotypes, stats, and set-extension entries.
4. Refpkg document: logical components, `refpkg_id`, validation failures, and why
   archive-byte hashes are not the identity.
5. Phylotype-set document: manifest, `set_id`/`parent_set_id`, dataset namespace,
   threshold replace/append behavior, stable IDs, new-group creation, and repeated
   extension.
6. Full composition examples remain, but are secondary to standalone entry examples.

## Verify

- Every command matches its entry's `--help` and schema.
- Examples use generic v2 paths and parameterized threshold discovery.
- No docs describe phylotypes as part of taxonomy or imply unassigned SVs.

## Commit / PR

- Title: `docs: modular refpkg-centered and versioned phylotype-set workflows`
