params.help = false
params.manifest = null
params.legacy_redup = false
params.output = '.'

nextflow.enable.dsl=2

include { read_manifest } from '../../modules/manifest'
include { output_failed; preprocess_wf } from '../../modules/preprocess'
include { dada2_wf } from '../../modules/dada2'
include { ExportSvPlacementInputs } from '../../modules/sv_h5ad'

def helpMessage() {
    log.info """
    ─────────────────────────────────────
    maliampi / sv
    ─────────────────────────────────────
    Generate sequence variants from raw reads using DADA2.

    Usage:
      nextflow run jgolob/maliampi/subworkflows/local/sv.nf [options]

    Required:
      --manifest      CSV file listing samples (columns: specimen, R1, R2;
                        optional: batch, I1, I2)

    Options:
      --legacy_redup  Export legacy pplacer-compatible weight files (default: false)
      --output        Output directory (default: .)
      --help          Show this help message
    """.stripIndent()
}

workflow sv {
    take:
    manifest_file

    main:
    manifest = read_manifest(manifest_file)

    preprocess_wf(
        manifest.valid_paired_indexed,
        manifest.valid_paired,
        manifest.valid_unpaired,
    )

    dada2_wf(
        preprocess_wf.out.miseq_pe,
        preprocess_wf.out.miseq_se,
        preprocess_wf.out.pyro,
    )
    ExportSvPlacementInputs(dada2_wf.out.sv_h5ad, params.legacy_redup)

    failures_for_report = manifest.other.map { row -> [row.specimen, 'failed at manifest'] }
        .mix(preprocess_wf.out.empty.map { row -> [row[0], 'preprocessing'] })
        .mix(dada2_wf.out.failures)
        .toList()
        .map { rows -> [
            rows.collect { row -> row[0] },
            rows.collect { row -> row[1] },
        ] }
    output_failed(failures_for_report)

    emit:
    sv_h5ad = dada2_wf.out.sv_h5ad
    sv_fasta = ExportSvPlacementInputs.out.sv_fasta
    sv_multiplicity = ExportSvPlacementInputs.out.multiplicity
    sv_weights = ExportSvPlacementInputs.out.sv_weights
    failures = dada2_wf.out.failures
}

workflow {
    if (params.help || params.manifest == null) {
        helpMessage()
        if (!params.help) {
            error "Missing required parameter(s): --manifest"
        }
        return
    }
    sv(channel.fromPath(params.manifest, checkIfExists: true))
}
