# Phase 21 — Pinned phylotypes v2.1.0 container

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/21-phylotypes-image` ·
**Files:** `docker/phylotypes/Dockerfile` · **Depends on:** 20a released

## Goal

Build a reproducible image containing both `phylotypes` and
`extend_phylotypes` from the reviewed v2.1.0 release.

## Steps

1. Pin installation to the exact v2.1.0 commit recorded in Phase 20a; do not use a
   moving branch or unverified local dirty checkout.
2. Build `golob/phylotypes:2.1.0` from a Python version supported by the upstream
   `pyproject.toml` and record the source commit as an OCI label.
3. Smoke-test:
   - `phylotypes --help`
   - `extend_phylotypes --help`
   - Initial mapping from fixture jplace
   - Extension in which one SV joins an old group and another creates a new group
4. Building locally and publishing are separate actions. Do not push the image
   without explicit approval.

## Verify

- The image reports phylotypes 2.1.0 and the pinned commit.
- Both CLIs produce complete, schema-valid mappings.
- Rebuilding from the same source produces functionally identical results.

## Commit / PR

- Title: `feat(container): pin phylotypes 2.1.0 generate/extend tools`
