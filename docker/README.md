# Imported and legacy container inventory

Every `golob/*` image is intentionally tracked here until it is replaced with a
current public BioContainer or rebuilt from imported source.  These images are
not considered current merely because their tag is pinned.

| Image | Replacement status | Source / action |
|---|---|---|
| `golob/fastatools:0.8.0A` | No BioContainers replacement found | Dockerfile imported from `jgolob/fastatools` commit `d049807`. Import its build-context scripts before rebuilding. |
| `golob/dada2-pplacer:0.8.0__bcw_0.3.1A` | No BioContainers replacement found | Dockerfile imported from `jgolob/dada2-pplacer` commit `c30d045`. Import its Python package before rebuilding; this is the DADA2 1.38 bridge risk. |
| `golob/barcodecop:0.5__bc_1` | No BioContainers replacement found | Source found at `jgolob/barcodecop` commit `61c468d`; it has no Dockerfile. Add a maintained Python-image recipe. |
| `golob/dada2-fast-combineseqtab:0.5.0__1.12.0__BCW_0.3.1` | No source located | Flag for source recovery or replacement before any rebuild. |
| `golob/goodsfilter:0.1.6` | No source located | Flag for source recovery or replacement before any rebuild. |
| `golob/seqinfo_taxonomy_sync:0.3.0` | No source located | Flag for source recovery or replacement before any rebuild. |
| `ghcr.io/jgolob/phylotypes:2.1.0` | Independent phylotypes utility | Pinned external image; do not vendor it into MaliAmPi. |

`golob/gappa:0.3`, `golob/pplacer:*`, and `golob/taxtastic:*` have been replaced
by current official BioContainers in `conf/base.config`.

## `maliampi-tools`

`ghcr.io/jgolob/maliampi-tools` is intentionally published for `linux/amd64`
only. `rds2py` provides the required Linux x86_64 wheel but no functioning
Linux ARM build. The `maliampi_tools` process label requests
`--platform linux/amd64`; Docker Desktop on Apple Silicon therefore runs this
one helper image through emulation. Build it explicitly with:

```sh
docker build --platform linux/amd64 -f docker/maliampi-tools/Dockerfile \
  -t ghcr.io/jgolob/maliampi-tools:0.1.0 .
```
o
