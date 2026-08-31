params.help = false
params.dedup_jplace = null
params.specimen_jplaces = null
params.output = '.'

nextflow.enable.dsl=2

include { GappaEDPL; GappaKRD; GappaEPCA } from '../../modules/stats'

def helpMessage() {
    log.info """
    ─────────────────────────────────────
    maliampi / stats
    ─────────────────────────────────────
    Compute placement statistics: EDPL, KR distance, and edge PCA.

    Usage:
      nextflow run jgolob/maliampi/subworkflows/local/stats.nf [options]

    Required:
      --dedup_jplace      Deduplicated jplace file from placement
      --specimen_jplaces  Per-specimen jplace files (glob or directory path ending in /)

    Options:
      --output            Output directory (default: .)
      --help              Show this help message
    """.stripIndent()
}

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

workflow {
    if (params.help || params.dedup_jplace == null || params.specimen_jplaces == null) {
        helpMessage()
        if (!params.help) {
            def missing = []
            if (params.dedup_jplace == null) missing << '--dedup_jplace'
            if (params.specimen_jplaces == null) missing << '--specimen_jplaces'
            error "Missing required parameter(s): ${missing.join(', ')}"
        }
        return
    }
    def specimen_glob = params.specimen_jplaces.endsWith('/') ? "${params.specimen_jplaces}*.jplace*" : params.specimen_jplaces
    stats(
        channel.fromPath(params.dedup_jplace, checkIfExists: true),
        channel.fromPath(specimen_glob, checkIfExists: true),
    )
}
