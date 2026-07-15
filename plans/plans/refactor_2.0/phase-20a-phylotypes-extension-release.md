# Phase 20a — Upstream phylotypes extension API and v2.1.0 release

**Flag:** 🔴 MAJOR · **Repository:** `jgolob/phylotypes` ·
**Branch:** `feature/extend-existing-set` · **Depends on:** none in MaliAmPi

## Goal

Provide a supported API/CLI for extending an existing threshold-specific phylotype
mapping. Every new SV joins exactly one existing or newly created phylotype; none
are orphaned, unassigned, or omitted.

## Required CLI/API

Add `extend_phylotypes` with inputs:

- Previous cumulative jplace
- Previous complete `phylotype,sv` mapping
- New jplace on the same reference tree
- Distance metric, distance threshold, LWR overlap, and distal-length policy
- Output path for the complete updated mapping

## Algorithm

1. Validate identical/canonically equivalent jplace trees and required fields.
2. Initialize the assignment pool from the existing mapping while retaining every
   existing phylotype ID and membership.
3. Assign compatible new SVs to existing groups using the configured metric and
   threshold.
4. Cluster all remaining new SVs among themselves at that same threshold and add
   new groups.
5. Allocate new IDs after the largest existing `pt__NNNNN`, deterministically.
6. Never merge, renumber, or delete an existing phylotype during extension.
7. Require the updated mapping to partition old+new SVs exactly once.

Also make initial `to_long()` ID allocation deterministic: sort groups by descending
size, then by their lexicographically smallest SV (and full sorted membership as the
final tie-break). A versioned set cannot depend on incidental cluster/list order.

Do not reuse the current `add_phylotypes.py` behavior that merely logs and omits
orphans. Reuse tested distance/pool primitives from the v2 incremental code, but do
not run its final reconciliation across established groups.

## Tests and release

- Existing IDs remain byte-for-byte stable.
- New SVs that fit join expected groups.
- Nonmatching new SVs form deterministic new phylotypes.
- Mixed existing/new inputs are a complete partition with no duplicates.
- Tree/field incompatibility fails before assignment.

Release as **v2.1.0**, record the immutable tag commit, and only then update the
MaliAmPi container pin. Publishing the upstream release requires explicit approval.

## PR

- Title: `feat: extend existing phylotype sets without unassigned SVs`
