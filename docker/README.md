# Container policy

MaliAmPi uses public BioContainers for external engines and one internal helper
image, `ghcr.io/jgolob/maliampi-tools`. Barcodecop is installed from PyPI into
that helper image. Phylotypes is independent and pinned separately in
`conf/base.config`.

The helper image is intentionally `linux/amd64` because the required `rds2py`
wheel is not functional on ARM. Build it explicitly with:

```sh
docker build --platform linux/amd64 -f docker/maliampi-tools/Dockerfile \
  -t ghcr.io/jgolob/maliampi-tools:0.1.0 .
```

No retired helper repository or legacy custom image is a source of truth.

## CI and releases

`.github/workflows/container.yaml` builds the AMD64 helper image on pull
requests and changes to `master`, but does not publish it. It imports the image
locally and verifies the installed H5AD, RDS, Good's, and Swarm command
interfaces. Its Dockerfile `test` target also runs the full Python test suite
inside the AMD64 build before the runtime image is accepted.

A `vX.Y.Z` tag runs the same validation and then publishes
`ghcr.io/jgolob/maliampi-tools:X.Y.Z`, a content-addressed `sha-<commit>` tag,
and the floating `latest` tag. Releases include build provenance and an SBOM,
and the tag also creates a GitHub Release `vX.Y.Z` with auto-generated notes
linking the published image tags. A `workflow_dispatch` run with `publish: true`
and an explicit `version` is available as a manual escape hatch.
