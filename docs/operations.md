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
relevant pull requests and on `master` changes without publishing. Its Docker
`test` target runs the Python suite inside the image; the runtime image then
receives CLI smoke checks.

A `vX.Y.Z` tag runs the same build-and-validate job and then publishes. It
pushes `ghcr.io/jgolob/maliampi-tools:X.Y.Z`, `sha-<commit>`, and the floating
`latest` tag — including SBOM and provenance attestations — and creates a GitHub
Release `vX.Y.Z` with auto-generated notes that link the published image tags.
Publication requires no manual gating; a `workflow_dispatch` run with
`publish: true` and an explicit `version` is available as a manual escape hatch.

## Verified smoke path

The bundled minimal data has been exercised through the helper image, standalone
Good's, standalone Swarm, DADA2 SV inference, refpkg construction, EPA-ng,
taxonomy, phylotypes, and stats. See `maliampi-practice-data/` for the fixture
inputs. Docker/Nextflow execution remains an explicit operational action; CI
validation on pull requests and `master` does not publish images. Only a
`vX.Y.Z` tag publishes and cuts a release.
