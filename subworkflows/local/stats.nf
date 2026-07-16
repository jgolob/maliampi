nextflow.enable.dsl=2

include { GappaEDPL; GappaKRD; GappaEPCA } from '../../modules/stats'

workflow stats {
    take:
    dedup_jplace
    specimen_jplaces

    main:
    GappaEDPL(dedup_jplace)
    specimen_jplaces.collect().set { collected_specimen_jplaces }
    GappaKRD(collected_specimen_jplaces)
    GappaEPCA(collected_specimen_jplaces)

    emit:
    edpl = GappaEDPL.out.edpl
    krd = GappaKRD.out.krd
    epca_projection = GappaEPCA.out.projection
    epca_transformation = GappaEPCA.out.transformation
}

workflow stats_entry {
    if (params.dedup_jplace == null || params.specimen_jplaces == null) {
        error 'stats requires --dedup_jplace and --specimen_jplaces (a glob or directory)'
    }
    def specimen_glob = params.specimen_jplaces.endsWith('/') ? "${params.specimen_jplaces}*.jplace*" : params.specimen_jplaces
    stats(
        channel.fromPath(params.dedup_jplace, checkIfExists: true),
        channel.fromPath(specimen_glob, checkIfExists: true),
    )
}
