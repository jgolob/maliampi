nextflow.enable.dsl=2

params.help = false
params.mode = null
params.sv_fasta = null
params.sv_long = null
params.sharetable = null
params.sv_map = null
params.sv_weights = null
params.dataset_id = null
params.dedup_jplace = null
params.output = '.'

include { NormalizeSvLong; SharetableToSvLong; MapWeightsToSvLong; ValidateSvMapWeights; PplacerReduplicate } from '../../modules/conversion'

def helpMessage() {
    log.info """
    ─────────────────────────────────────
    maliampi / convert
    ─────────────────────────────────────
    Convert legacy SV abundance formats into canonical H5AD, or
    reduplicate a dedup jplace file.

    Usage:
      nextflow run jgolob/maliampi/subworkflows/local/convert.nf [options]

    Mode 1 — Normalize legacy SV data to H5AD (default when --sv_fasta given):
      Supply --sv_fasta and exactly one abundance source:

        --sv_fasta + --sv_long                    Long-format CSV (specimen,sv,count)
        --sv_fasta + --sharetable                 Mothur-style shared table
        --sv_fasta + --sv_map + --sv_weights      Legacy map + weights pair
        --sv_fasta + --sv_long + --sv_map + --sv_weights   All four (cross-validates)

      Optional:
        --dataset_id    Dataset identifier (default: empty)

    Mode 2 — Pplacer reduplication (when --dedup_jplace given):
      --dedup_jplace    Deduplicated jplace file
      --sv_weights      SV weights file for reduplication

    Common options:
      --output          Output directory (default: .)
      --help            Show this help message

    Tip: For quick one-off conversion without Nextflow:
      maliampi-sv-import --fasta seqs.fasta --long counts.csv \\
        --project-id myproject --dataset-id mydataset --output-h5ad sv.h5ad
    """.stripIndent()
}

workflow normalize_sv {
    take:
    sv_fasta
    sv_long
    dataset_id

    main:
    NormalizeSvLong(sv_fasta, sv_long, dataset_id)

    emit:
    sv_h5ad = NormalizeSvLong.out.sv_h5ad
}

workflow {
    if (params.help) {
        helpMessage()
        return
    }

    // Auto-detect mode from supplied params
    def is_normalize = params.sv_fasta != null
    def is_redup = params.dedup_jplace != null && params.sv_fasta == null

    if (!is_normalize && !is_redup) {
        helpMessage()
        error "Supply --sv_fasta (normalize mode) or --dedup_jplace (reduplication mode). Run with --help for details."
    }

    if (is_redup) {
        // Pplacer reduplication mode
        if (params.sv_weights == null) {
            helpMessage()
            error "Missing required parameter: --sv_weights (required for reduplication mode)"
        }
        PplacerReduplicate(
            channel.fromPath(params.dedup_jplace, checkIfExists: true),
            channel.fromPath(params.sv_weights, checkIfExists: true),
        )
    } else {
        // Normalize mode — detect abundance source
        def modes = [
            params.sharetable != null,
            params.sv_map != null && params.sv_weights != null && params.sv_long == null,
            params.sv_long != null && params.sv_map == null && params.sv_weights == null,
            params.sv_long != null && params.sv_map != null && params.sv_weights != null,
        ]
        if (modes.count { it } != 1) {
            helpMessage()
            error 'Supply exactly one legacy abundance source with --sv_fasta. Run with --help for details.'
        }
        def dataset_id = params.dataset_id ?: ''

        if (params.sharetable != null) {
            SharetableToSvLong(channel.fromPath(params.sharetable, checkIfExists: true))
            normalize_sv(channel.fromPath(params.sv_fasta, checkIfExists: true), SharetableToSvLong.out.long, dataset_id)
        } else if (params.sv_map != null && params.sv_weights != null && params.sv_long == null) {
            MapWeightsToSvLong(
                channel.fromPath(params.sv_map, checkIfExists: true),
                channel.fromPath(params.sv_weights, checkIfExists: true),
            )
            normalize_sv(channel.fromPath(params.sv_fasta, checkIfExists: true), MapWeightsToSvLong.out.long, dataset_id)
        } else if (params.sv_map != null && params.sv_weights != null) {
            ValidateSvMapWeights(
                channel.fromPath(params.sv_fasta, checkIfExists: true),
                channel.fromPath(params.sv_long, checkIfExists: true),
                channel.fromPath(params.sv_map, checkIfExists: true),
                channel.fromPath(params.sv_weights, checkIfExists: true),
            )
            normalize_sv(ValidateSvMapWeights.out.fasta, ValidateSvMapWeights.out.long, dataset_id)
        } else {
            normalize_sv(
                channel.fromPath(params.sv_fasta, checkIfExists: true),
                channel.fromPath(params.sv_long, checkIfExists: true),
                dataset_id,
            )
        }
    }
}
