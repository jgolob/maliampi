# Phase 09a — Refpkg identity and validation

**Flag:** 🔴 MAJOR · **Branch:** `refactor2/09a-refpkg-validation` ·
**Files:** `modules/refpkg_validate.nf`, `subworkflows/local/validate_refpkg.nf` ·
**Depends on:** 09

## Goal

Make the refpkg the immutable compatibility anchor. Every workflow that consumes an
existing refpkg must validate it and carry a content-derived `refpkg_id` before any
expensive downstream work starts.

## Identity contract

Compute `refpkg_id` as SHA-256 over sorted logical component names and each
component's SHA-256. Do not hash tar ordering, timestamps, compression bytes, logs,
or rollback metadata. Repacking identical logical content must preserve the ID.

## Validation

`ValidateRefpkg` must:

- Reject unsafe archive members, duplicate basenames, missing components, and
  unsupported format/locus metadata.
- Verify legacy `CONTENTS.json` MD5 entries where present, then compute SHA-256 for
  every logical component.
- Require matching FASTA/Stockholm reference identifiers and alignment lengths.
- Require tree leaves to equal alignment identifiers.
- Require sequence-info and taxonomy coverage for all tree leaves.
- Require nonempty, parseable model and covariance-model components.
- Emit `refpkg.identity.json`, `refpkg.validation.json`, and the original refpkg.

Expose imported workflow `validate_refpkg(refpkg)` and standalone
`validate_refpkg_entry --refpkg`. The imported workflow emits a validated tuple
that downstream stages consume; do not pass an unvalidated archive onward.

## Fail-fast gating

For workflows starting from an existing refpkg, derive the manifest/input channel
from the validation output token so FASTQ preprocessing cannot begin concurrently
with a validation that may fail.

## Verify

- The fixture refpkg validates and receives a stable ID across equivalent repacks.
- Tampering with every component type fails checksum validation.
- Tree/alignment/taxonomy mismatches fail with a specific diagnostic.
- Unsafe tar members never extract outside the task directory.

## Commit / PR

- Title: `feat(refpkg): content identity and fail-fast validation entry`
