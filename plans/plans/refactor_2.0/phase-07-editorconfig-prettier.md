# Phase 07 — `.editorconfig` + Prettier

**Flag:** 🟢 SIMPLE · **Branch:** `refactor2/07-editorconfig-prettier` ·
**Files:** `.editorconfig`, `.prettierrc.yml`, `.prettierignore` (all new) ·
**Depends on:** none (can run any time)

## Goal
Add the nf-core-standard formatting baseline so style stops drifting.

## Steps
1. **`.editorconfig`** (nf-core standard):
   ```editorconfig
   root = true

   [*]
   charset = utf-8
   end_of_line = lf
   insert_final_newline = true
   trim_trailing_whitespace = true
   indent_size = 4
   indent_style = space

   [*.{yml,yaml}]
   indent_size = 2

   [*.{md,nf,config}]
   trim_trailing_whitespace = false   # markdown line breaks; nf heredocs

   [*.{json}]
   insert_final_newline = unset
   ```
2. **`.prettierrc.yml`**:
   ```yaml
   printWidth: 120
   tabWidth: 4
   ```
3. **`.prettierignore`** — exclude data/fixtures and generated dirs:
   ```
   tests/
   maliampi-practice-data/
   data/
   work/
   output/
   *.nf            # nextflow not handled by prettier
   ```

## Verify
```bash
npx prettier --check "**/*.{yml,yaml,json,md}" || true   # lists files; fix obvious ones
```
Do **not** mass-reformat existing files in this PR — just add the config. A later
optional cleanup can run `prettier --write`.

## Commit / PR
- Title: `chore: add editorconfig + prettier config`
