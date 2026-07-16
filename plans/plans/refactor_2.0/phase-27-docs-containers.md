# Phase 27 — `docs/containers.md`

**Flag:** 🟢 SIMPLE · **Branch:** `refactor2/27-docs-containers` ·
**Files:** `docs/containers.md` (new) · **Depends on:** 21, 23 (for final pins)

## Goal
A single tracked table of every container the pipeline uses: image, pin, whether a
biocontainer equivalent exists, and how to rebuild bespoke ones.

## Content
1. Table columns: **tool · image:tag · source (biocontainer/golob) · equivalent?
   · rebuild Dockerfile**.
2. Rows for every `params.container__*` across `modules/` (grep them):
   ```bash
   git grep -hoE 'params\.container__[a-zA-Z0-9_]+ *= *['"'"'"][^'"'"'"]+['"'"'"]' modules/ | sort -u
   ```
3. Bespoke `golob/*` images with **no** biocontainer equivalent — note their purpose
   and link a real reproducible build source where one exists. Do not create empty
   or placeholder Dockerfiles merely to fill the table: `fastatools`,
   `dada2-pplacer`, `dada2-fast-combineseqtab`,
   `goodsfilter`, `seqinfo_taxonomy_sync`, `pplacer`, `phylotypes`.
4. **Follow-up section**: images still worth migrating to biocontainers later
   (`taxtastic`, `barcodecop`) — note candidate biocontainer + open question.

## Verify
- Every `container__*` in the codebase appears as a row (diff the grep list vs the
  table).

## Commit / PR
- Title: `docs: container inventory + provenance (docs/containers.md)`
