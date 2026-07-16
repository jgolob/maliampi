nextflow.enable.dsl=2

include { NormalizeSvLong; SharetableToSvLong; MapWeightsToSvLong; ValidateSvMapWeights; PplacerReduplicate } from '../../modules/conversion'

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

workflow normalize_sv_entry {
    def modes = [
        params.sv_fasta != null && params.sharetable != null,
        params.sv_fasta != null && params.sv_map != null && params.sv_weights != null && params.sv_long == null,
        params.sv_fasta != null && params.sv_long != null && params.sv_map == null && params.sv_weights == null,
        params.sv_fasta != null && params.sv_long != null && params.sv_map != null && params.sv_weights != null,
    ]
    if (modes.count { it } != 1) {
        error 'Supply exactly one legacy import source: --sv_fasta --sharetable; --sv_fasta --sv_map --sv_weights; --sv_fasta --sv_long; or all four validation files'
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

workflow pplacer_reduplicate_entry {
    if (params.dedup_jplace == null || params.sv_weights == null) {
        error 'pplacer reduplication requires --dedup_jplace and --sv_weights'
    }
    PplacerReduplicate(
        channel.fromPath(params.dedup_jplace, checkIfExists: true),
        channel.fromPath(params.sv_weights, checkIfExists: true),
    )
}
