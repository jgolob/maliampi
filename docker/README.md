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
requests, changes to `master`, and `pre-v*` milestone tags, but does not publish it. It imports the image
locally and verifies the installed H5AD, RDS, Good's, and Swarm command
interfaces. Its Dockerfile `test` target also runs the full Python test suite
inside the AMD64 build before the runtime image is accepted.

Publishing to GHCR is deliberately disabled until a repository administrator:

1. creates the `container-release` GitHub Environment;
2. restricts it to protected `tools-v*` tags and requires reviewers; and
3. sets repository variable `CONTAINER_RELEASE_ENABLED` to `true`.

After that, a protected `tools-vX.Y.Z` tag or an approved manual dispatch can
publish `ghcr.io/jgolob/maliampi-tools:X.Y.Z` and a content-addressed
`sha-<commit>` tag. Releases include build provenance and an SBOM. There is no
automatic `latest` or `edge` tag. A `pre-vX.Y.Z` tag is a validation milestone
only; it cannot publish an image.
