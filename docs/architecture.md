# Architecture

MaliAmPi is a refpkg-centered, modular amplicon workflow. The composed workflow
is convenient, but every major stage can consume and produce durable artifacts
without rerunning raw-read processing.

```text
FASTQ ──> SV inference ──> sv.h5ad ──> placement ──> taxonomy H5AD
                              │              ├──> phylotype H5AD
                              │              └──> placement statistics
                              └──> refpkg registry / refpkg construction
```

## Canonical SV contract

`sv.h5ad` is the only internal specimen-by-SV abundance contract.

| Location | Contents |
|---|---|
| `X` | Integer, nonnegative specimen × SV CSR counts; zeros are omitted. |
| `obs` | Exact `project_id` and deposited specimen ID, a full identity SHA-256, and run/batch provenance. |
| `var` | Composite SV ID, uppercase IUPAC sequence, full sequence SHA-256, length, and total abundance. |
| `uns["maliampi"]` | Schema version, project/dataset/tool provenance, and optional refpkg identity. |

SV sequence identity is the full SHA-256 of the normalized uppercase IUPAC DNA
sequence. The readable `sv-<decimal>-<16-hex>` identifier is derived from that
digest, but is not authoritative on its own. Observations aggregate repeated
runs of the same project/specimen into one row.

Every helper that reads an SV H5AD validates dimensions, count type,
non-negativity, identity uniqueness, hashes, composite IDs, and the absence of
empty observations or variables.

## Boundary formats

FASTA, gappa multiplicity CSV, and optional pplacer map/weights are generated
only at external-tool boundaries. Legacy FASTA/long/share/map/weights inputs
are supported by the conversion entrypoint, which immediately writes a
validated H5AD. They are not persistent internal workflow contracts.

Taxonomy and SV-to-phylotype mappings are Parquet. Taxonomy-rank and
phylotype-threshold abundance artifacts are H5AD files with integer `X` and a
sparse `relative_abundance` layer.

## Refpkg integration

The refpkg includes an SV identity registry. Before placement, MaliAmPi checks
the input H5AD against that registry and records immutable SV-artifact and
refpkg identities in downstream provenance. Existing sequence hashes resolve
to known SVs; genuinely new sequences are placed and may extend the phylotype
set. A mismatched SV identity is rejected rather than silently remapped.
