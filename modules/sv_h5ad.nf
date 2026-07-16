nextflow.enable.dsl=2

/* Canonical SV artifact adapters. CSV output is deliberately boundary-only. */

process FinalizeSvH5ad {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
    path pre_chimera_h5ad
    path retained_sequences

    output:
    path 'sv.h5ad', emit: sv_h5ad

    script:
    """
    maliampi-sv-filter --sequences ${pre_chimera_h5ad} ${retained_sequences} sv.h5ad
    """
}

process ExportSvPlacementInputs {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path sv_h5ad
    val legacy_redup

    output:
    path 'placement/sv.fasta', emit: sv_fasta
    path 'placement/sv.multiplicity.csv', emit: multiplicity
    path 'placement/sv.weights.csv', optional: true, emit: sv_weights

    script:
    def redup = legacy_redup ? '--pplacer-reduplication' : ''
    """
    maliampi-sv-export ${sv_h5ad} placement ${redup}
    """
}

process WriteSvRegistry {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path sv_h5ad

    output:
    path 'sv_registry.parquet', emit: registry

    script:
    """
    maliampi-refpkg-registry ${sv_h5ad} sv_registry.parquet
    """
}

process ValidateSvRegistry {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path sv_h5ad
    path sv_registry

    output:
    path 'sv_registry.validation.json', emit: validation

    script:
    """
    maliampi-refpkg-validate ${sv_h5ad} ${sv_registry} > sv_registry.validation.json
    """
}
