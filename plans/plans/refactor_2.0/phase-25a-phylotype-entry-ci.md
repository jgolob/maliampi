# Phase 25a — Phylotypes and set-extension entry CI

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/25a-phylotype-ci` ·
**Files:** `.github/workflows/test.yaml`, phylotype fixtures/assertions ·
**Depends on:** 22a, 24

## Goal

Add CI for the parallel phylotypes pathway and repeatable extension of a validated,
versioned set without coupling this work to monolith deletion.

## Jobs

1. `phylotypes_entry`: parameterized threshold list, mappings, count/RA tables, and
   initial set manifest.
2. Threshold semantics: replace defaults, append values, sort/deduplicate, and reject
   invalid/KR-without-calibration input.
3. `extend_phylotype_set_entry`: one SV joins an old group, one creates a new group,
   all existing IDs remain stable, and every new SV appears once.
4. New-threshold extension: full cumulative mapping is created while stored
   unrequested thresholds remain present.
5. Repeated extension: new `set_id`, correct `parent_set_id`, cumulative jplace, and
   count conservation.
6. Negative cases: refpkg/tree mismatch, tampered bundle, duplicate dataset/SV IDs.
7. Full composition default/skip behavior after phylotypes is added to main.

Discover outputs from the effective threshold list; never hardcode three directories.

## Commit / PR

- Title: `test(ci): multi-threshold phylotypes and versioned-set extension`
