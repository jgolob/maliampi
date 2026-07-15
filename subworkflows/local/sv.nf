nextflow.enable.dsl=2

include { read_manifest } from '../../modules/manifest'
include { output_failed; preprocess_wf } from '../../modules/preprocess'
include { dada2_wf } from '../../modules/dada2'
include { ExportSvPlacementInputs } from '../../modules/sv_h5ad'

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

workflow sv_entry {
    if (params.manifest == null) {
        error 'Missing required --manifest'
    }
    sv(channel.fromPath(params.manifest, checkIfExists: true))
}
