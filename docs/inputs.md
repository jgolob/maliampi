# Inputs and identities

## Manifest

The FASTQ manifest is CSV with these columns:

| Column | Required | Meaning |
|---|---:|---|
| `specimen` | Yes | Deposited specimen identifier, preserved exactly. |
| `batch` | Recommended | Extraction/sequencing run provenance. |
| `R1` | Yes | Forward FASTQ path. |
| `R2` | Paired reads | Reverse FASTQ path. |
| `I1`, `I2` | Indexed reads | Index FASTQs for barcodecop verification. |
| `data_type` | Optional | Set to `pyro` or `16S_pyro` for pyrosequencing input. |

Use paths that are valid from the directory where Nextflow is launched. A row
that cannot be read is reported in `failed_specimens.csv`; it is not silently
incorporated into the H5AD.

## Project and dataset IDs

`--project_id` and `--dataset_id` are both required for active SV workflows.

- `project_id` identifies the deposited study/project and participates in each
  observation identity.
- `dataset_id` identifies this processed SV artifact and is recorded in H5AD
  provenance and phylotype-set metadata.

Do not reuse a project/specimen pair for unrelated biological material. Repeated
technical runs or batches of the same deposited specimen intentionally
aggregate into one H5AD observation, with all contributing batches retained in
provenance.

## Refpkg inputs

Build a refpkg from:

- `--repo_fasta`: reference sequences;
- `--repo_si`: reference sequence-information table;
- `--email`: contact address required for the refpkg workflow; and optionally
- `--taxdmp`: a local NCBI taxonomy dump.

For routine new FASTQ data, prefer `--refpkg existing-refpkg.tar.gz`. This
keeps placement, taxonomy, and phylotypes anchored to the existing reference
package rather than rebuilding it for every run.

## Legacy SV imports

The conversion entrypoint accepts exactly one of the following source sets:

1. `--sv_fasta` and `--sharetable`;
2. `--sv_fasta`, `--sv_map`, and `--sv_weights`;
3. `--sv_fasta` and `--sv_long`; or
4. all of FASTA, long, map, and weights, which are cross-validated before
   import.

Each import immediately produces validated `sv.h5ad` rather than continuing as
a CSV workflow.
