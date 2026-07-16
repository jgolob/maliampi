# Operations, containers, and releases

## Resource configuration

Resource assignments belong in Nextflow configuration labels, not in workflow
processes. The shipped labels are `io_limited`, `io_net`, `multithread`,
`mem_medium`, and `mem_veryhigh`. Copy and adjust the profile for the target
executor rather than editing modules to change memory or CPU allocations.

`conf/test.config` is intentionally small and exists for the bundled smoke
dataset. It is not a production resource profile.

## Containers

External engines use pinned BioContainers. `maliampi-tools` is the unified
helper image; barcodecop is installed from PyPI in that image. Phylotypes is an
independent, separately released image.

The helper is intentionally AMD64-only because `rds2py` lacks a functioning ARM
build. The `maliampi_tools` process label requests `--platform linux/amd64`, so
Docker Desktop on Apple Silicon uses emulation.

Build and test the helper image locally:

```sh
docker build --platform linux/amd64 --target test \
  -f docker/maliampi-tools/Dockerfile .
docker build --platform linux/amd64 \
  -t maliampi-tools:local -f docker/maliampi-tools/Dockerfile .
```

## CI and GHCR

`.github/workflows/container.yaml` builds and tests the AMD64 helper image on
relevant pull requests, `master` changes, and `pre-v*` milestone tags without
publishing. Its Docker
`test` target runs the Python suite inside the image; the runtime image then
receives CLI smoke checks.

GHCR publication is disabled by default. To enable approved releases, a
repository administrator must create the `container-release` environment,
restrict it to protected `tools-v*` tags, require reviewers, and set repository
variable `CONTAINER_RELEASE_ENABLED=true`.

An approved `tools-vX.Y.Z` tag publishes
`ghcr.io/jgolob/maliampi-tools:X.Y.Z` and `sha-<commit>`, including SBOM and
provenance attestations. A `pre-vX.Y.Z` tag is validation-only and cannot
publish. No floating `latest` or `edge` tag is published.

## Verified smoke path

The bundled minimal data has been exercised through the helper image, standalone
Good's, standalone Swarm, DADA2 SV inference, refpkg construction, EPA-ng,
taxonomy, phylotypes, and stats. See `maliampi-practice-data/` for the fixture
inputs. Docker/Nextflow execution remains an explicit operational action; CI
validation does not publish images or run a remote release.
