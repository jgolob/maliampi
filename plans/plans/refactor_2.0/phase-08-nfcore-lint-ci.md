# Phase 08 — `.nf-core.yml` + CI `lint` job

**Flag:** 🟡 MEDIUM · **Branch:** `refactor2/08-lint-ci` ·
**Files:** `.nf-core.yml` (new), `.github/workflows/lint.yml` (new) ·
**Depends on:** 07

## Goal
Run automated linting in CI — **non-blocking initially** (report-only), so we can
triage the findings into a punch-list without blocking merges.

## Steps
1. **`.nf-core.yml`** — declare this is a (non-canonical) nf-core-style pipeline and
   disable the checks that don't apply yet, so the linter is usable. Start minimal:
   ```yaml
   repository_type: pipeline
   lint:
     # Disable checks we haven't satisfied yet; re-enable as phases land.
     files_exist: false
     nextflow_config: false
     files_unchanged: false
   ```
   (Tune the disabled list after the first run; the goal is a green-ish baseline.)
2. **`.github/workflows/lint.yml`** — a separate workflow:
   ```yaml
   name: lint
   on: [push, pull_request, workflow_dispatch]
   jobs:
     editorconfig:
       runs-on: ubuntu-latest
       steps:
         - uses: actions/checkout@v4
         - uses: editorconfig-checker/action-editorconfig-checker@main
         - run: editorconfig-checker
     prettier:
       runs-on: ubuntu-latest
       steps:
         - uses: actions/checkout@v4
         - uses: actions/setup-node@v4
         - run: npx prettier --check "**/*.{yml,yaml,json,md}"
     nf-core:
       runs-on: ubuntu-latest
       continue-on-error: true        # non-blocking for now
       steps:
         - uses: actions/checkout@v4
         - uses: actions/setup-python@v5
         - run: pip install nf-core
         - run: nf-core pipelines lint || true
   ```

## Verify
- Push the branch; confirm the `lint` workflow appears and runs. `editorconfig` and
  `prettier` should pass (Phase 07 added the configs); `nf-core` may be red but is
  `continue-on-error`.

## Commit / PR
- Title: `ci: add non-blocking lint workflow (editorconfig, prettier, nf-core)`
- In the PR body, paste the nf-core lint summary as the **triage punch-list**.
