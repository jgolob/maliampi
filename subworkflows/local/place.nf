nextflow.enable.dsl=2

params.help = false
params.sv_h5ad = null
params.refpkg = null
params.placer = 'epang'
params.legacy_redup = false
params.epang_chunk_size = 5000
params.output = '.'

include { ValidateRefpkg } from '../../modules/refpkg_validate'
include { ExtractRefpkg } from '../../modules/refpkg_utils'
include { AlignSV; CombineAln_SV_refpkg; ConvertAlnToFasta; EPAngSplit; EPAngPlaceChunk; MergeJplace; PplacerPlacement; GappaSplit; PlacementIdentity } from '../../modules/place'
include { PplacerReduplicate } from '../../modules/conversion'
include { ExportSvPlacementInputs; ValidateSvRegistry } from '../../modules/sv_h5ad'

def helpMessage() {
    log.info """
    ─────────────────────────────────────
    maliampi / place
    ─────────────────────────────────────
    Place sequence variants onto a reference phylogeny.

    Usage:
      nextflow run jgolob/maliampi/subworkflows/local/place.nf [options]

    Required:
      --sv_h5ad       H5AD file with sequence variants
      --refpkg        Reference package (.tar.gz)

    Options:
      --placer        Placement engine: epang (default) or pplacer
      --legacy_redup      Emit legacy pplacer-reduplicated jplace (default: false)
      --epang_chunk_size  Query sequences per EPA-ng chunk (default: 5000)
      --output            Output directory (default: .)
      --help              Show this help message
    """.stripIndent()
}

workflow place {
    take:
    sv_h5ad
    refpkg
    placer
    legacy_redup

    main:
    if (!(placer in ['epang', 'pplacer'])) {
        error "Unsupported placement engine: ${placer}; use epang or pplacer"
    }
    ValidateRefpkg(refpkg)
    ExtractRefpkg(ValidateRefpkg.out.refpkg)
    ValidateSvRegistry(sv_h5ad, ExtractRefpkg.out.sv_registry)
    sv_h5ad.combine(ValidateSvRegistry.out.validation)
        .map { artifact, _validation -> artifact }
        .set { validated_sv_h5ad }
    ExportSvPlacementInputs(validated_sv_h5ad, legacy_redup)
    AlignSV(ExportSvPlacementInputs.out.sv_fasta, ExtractRefpkg.out.cm)
    CombineAln_SV_refpkg(AlignSV.out.stockholm, ExtractRefpkg.out.ref_aln_sto)

    if (placer == 'epang') {
        ConvertAlnToFasta(CombineAln_SV_refpkg.out.stockholm)
        EPAngSplit(ExtractRefpkg.out.ref_aln_fasta, ConvertAlnToFasta.out.fasta)
        // Scatter query sequences into chunks for parallel placement
        query_chunks = EPAngSplit.out.query.splitFasta(
            by: params.epang_chunk_size as int, file: true
        )
        EPAngSplit.out.reference
            .combine(query_chunks)
            .combine(ExtractRefpkg.out.model)
            .combine(ExtractRefpkg.out.tree)
            .set { placement_inputs }
        // placement_inputs is now [reference, chunk, model, tree] per chunk
        EPAngPlaceChunk(placement_inputs)
        MergeJplace(EPAngPlaceChunk.out.jplace.collect())
        dedup_jplace = MergeJplace.out.dedup_jplace
    } else {
        PplacerPlacement(CombineAln_SV_refpkg.out.stockholm, ValidateRefpkg.out.refpkg)
        dedup_jplace = PplacerPlacement.out.dedup_jplace
    }

    GappaSplit(dedup_jplace, ExportSvPlacementInputs.out.multiplicity)
    PlacementIdentity(dedup_jplace, ValidateRefpkg.out.identity, validated_sv_h5ad, placer)
    if (legacy_redup) {
        PplacerReduplicate(dedup_jplace, ExportSvPlacementInputs.out.sv_weights)
    }

    emit:
    dedup_jplace = dedup_jplace
    specimen_jplaces = GappaSplit.out.specimen_jplaces
    placement_identity = PlacementIdentity.out.identity
    validated_refpkg = ValidateRefpkg.out.refpkg
    refpkg_identity = ValidateRefpkg.out.identity
    redup_jplace = legacy_redup ? PplacerReduplicate.out.redup_jplace : channel.empty()
}

workflow {
    if (params.help || params.sv_h5ad == null || params.refpkg == null) {
        helpMessage()
        if (!params.help) {
            def missing = []
            if (params.sv_h5ad == null) missing << '--sv_h5ad'
            if (params.refpkg == null) missing << '--refpkg'
            error "Missing required parameter(s): ${missing.join(', ')}"
        }
        return
    }
    place(
        channel.fromPath(params.sv_h5ad, checkIfExists: true),
        channel.fromPath(params.refpkg, checkIfExists: true),
        params.placer,
        params.legacy_redup,
    )
}
