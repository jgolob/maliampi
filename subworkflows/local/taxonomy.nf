nextflow.enable.dsl=2

include { ValidateRefpkg } from '../../modules/refpkg_validate'
include { ExtractRefpkg } from '../../modules/refpkg_utils'
include { MakeGappaTaxonFile; GappaTaxonomy; ExtractGappaTaxonomy; AggregateTaxonomyRank } from '../../modules/taxonomy'

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

workflow taxonomy_entry {
    if (params.dedup_jplace == null || params.refpkg == null || params.sv_h5ad == null || params.placement_identity == null) {
        error 'taxonomy requires --dedup_jplace, --refpkg, --sv_h5ad, and --placement_identity'
    }
    taxonomy(
        channel.fromPath(params.dedup_jplace, checkIfExists: true),
        channel.fromPath(params.refpkg, checkIfExists: true),
        channel.fromPath(params.sv_h5ad, checkIfExists: true),
        channel.fromPath(params.placement_identity, checkIfExists: true),
    )
}
