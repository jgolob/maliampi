nextflow.enable.dsl=2

params.help = false
params.dedup_jplace = null
params.refpkg = null
params.sv_h5ad = null
params.placement_identity = null
params.output = '.'

include { ValidateRefpkg } from '../../modules/refpkg_validate'
include { ExtractRefpkg } from '../../modules/refpkg_utils'
include { MakeGappaTaxonFile; GappaTaxonomy; ExtractGappaTaxonomy; AggregateTaxonomyRank } from '../../modules/taxonomy'

def helpMessage() {
    log.info """
    ─────────────────────────────────────
    maliampi / taxonomy
    ─────────────────────────────────────
    Assign taxonomy to placed sequence variants using gappa.

    Usage:
      nextflow run jgolob/maliampi/subworkflows/local/taxonomy.nf [options]

    Required:
      --dedup_jplace        Deduplicated jplace file from placement
      --refpkg              Reference package (.tar.gz)
      --sv_h5ad             H5AD file with sequence variants
      --placement_identity  Placement identity JSON from placement step

    Options:
      --output              Output directory (default: .)
      --help                Show this help message
    """.stripIndent()
}

workflow taxonomy {
    take:
    dedup_jplace
    refpkg
    sv_h5ad
    placement_identity

    main:
    ValidateRefpkg(refpkg)
    ExtractRefpkg(ValidateRefpkg.out.refpkg)
    MakeGappaTaxonFile(ExtractRefpkg.out.leaf_info, ExtractRefpkg.out.taxonomy)
    GappaTaxonomy(dedup_jplace, MakeGappaTaxonFile.out.taxon_file)
    ExtractGappaTaxonomy(GappaTaxonomy.out.per_query, ExtractRefpkg.out.taxonomy, sv_h5ad)
    sv_h5ad.combine(ExtractGappaTaxonomy.out.sv_taxonomy)
        .combine(channel.of('phylum', 'class', 'order', 'family', 'genus', 'species'))
        .combine(ValidateRefpkg.out.identity)
        .combine(placement_identity)
        .set { table_inputs }
    AggregateTaxonomyRank(table_inputs)

    emit:
    sv_taxonomy = ExtractGappaTaxonomy.out.sv_taxonomy
    abundance_h5ad = AggregateTaxonomyRank.out.abundance_h5ad
    refpkg_identity = ValidateRefpkg.out.identity
}

workflow {
    def required = [
        '--dedup_jplace': params.dedup_jplace,
        '--refpkg': params.refpkg,
        '--sv_h5ad': params.sv_h5ad,
        '--placement_identity': params.placement_identity,
    ]
    def missing = required.findAll { _k, v -> v == null }.collect { k, _v -> k }

    if (params.help || missing) {
        helpMessage()
        if (!params.help) {
            error "Missing required parameter(s): ${missing.join(', ')}"
        }
        return
    }
    taxonomy(
        channel.fromPath(params.dedup_jplace, checkIfExists: true),
        channel.fromPath(params.refpkg, checkIfExists: true),
        channel.fromPath(params.sv_h5ad, checkIfExists: true),
        channel.fromPath(params.placement_identity, checkIfExists: true),
    )
}
