nextflow.enable.dsl=2

include { ValidateRefpkg } from '../../modules/refpkg_validate'
include { ExtractRefpkg } from '../../modules/refpkg_utils'
include { AlignSV; CombineAln_SV_refpkg; ConvertAlnToFasta; EPAngPlacement; PplacerPlacement; MakeSplit; GappaSplit; PlacementIdentity } from '../../modules/place'
include { PplacerReduplicate } from '../../modules/conversion'
include { ExportSvPlacementInputs; ValidateSvRegistry } from '../../modules/sv_h5ad'

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
        EPAngPlacement(
            ExtractRefpkg.out.ref_aln_fasta,
            ConvertAlnToFasta.out.fasta,
            ExtractRefpkg.out.model,
            ExtractRefpkg.out.tree,
        )
        dedup_jplace = EPAngPlacement.out.dedup_jplace
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

workflow place_entry {
    if (params.sv_h5ad == null || params.refpkg == null) {
        error 'place requires --sv_h5ad and --refpkg'
    }
    def placer = params.placer ?: 'epang'
    def legacy_redup = params.legacy_redup ?: false
    place(
        channel.fromPath(params.sv_h5ad, checkIfExists: true),
        channel.fromPath(params.refpkg, checkIfExists: true),
        placer,
        legacy_redup,
    )
}
