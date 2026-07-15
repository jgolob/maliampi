# Phase 23 — Biocontainer sweep + interval-KRD gappa migration

**Flag:** 🟡 MEDIUM (🔴 if outputs shift) · **Branch:** `refactor2/23-container-sweep`
· **Files:** all modules with `container__*`; `tests/` if fixtures change ·
**Depends on:** 20 (do DADA2 separately first)

## Goal
Bring every remaining biocontainer to current stable, retire the personal
`golob/gappa:0.3` in favor of the biocontainer, and pin everything to a build hash.

## Targets (verify each latest tag at implementation, then pin the hash)
| var / tool | current | →action |
|---|---|---|
| `container__infernal` | 1.1.4 / 1.1.5 split | unify on latest 1.1.5 build |
| `container__easel` | 0.47 | latest |
| `container__vsearch` | 2.22.1 | latest |
| trim-galore | 0.6.6 | latest |
| fastqc | v0.11.9 | latest 0.12.x |
| swarm | 3.1.2 | latest |
| `container__raxmlng` | 1.2.2 | confirm latest |
| raxml (og) | 8.2.4 | confirm latest 8.2.x |
| `container__epang` | 0.3.8 | confirm latest build |
| `container__gappa` (`interval_krd.nf`) | `golob/gappa:0.3` | → `quay.io/biocontainers/gappa:<latest>` |

Find each tag:
```bash
curl -s "https://quay.io/api/v1/repository/biocontainers/<tool>/tag/?limit=10&onlyActiveTags=true" \
  | grep -oE '"name": "[^"]*"'
```

## Steps
1. Bump one tool per commit (easier to bisect if a bump breaks output). Pin the full
   `tool:ver--buildhash`.
2. **interval-KRD gappa migration:** replace `golob/gappa:0.3` with the biocontainer
   already used elsewhere; confirm `interval_krd.nf`'s gappa invocation flags are
   compatible with the newer gappa (CLI changed between 0.3 and 0.9).
3. After each output-affecting bump, run the relevant workflow and refresh the
   matching `tests/` fixture **in the same commit**, explaining the change.

## Verify
```bash
git grep -n "golob/gappa" modules/      # empty after migration
nextflow run . -profile test,docker     # full auto run green on bumped images
```

## Commit / PR
- One PR, but **one commit per tool bump**. Title:
  `chore(containers): bump biocontainers to current stable + drop golob/gappa`
- 🔴 If any bump changes outputs, call it out per-commit with the fixture diff.
