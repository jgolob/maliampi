nextflow.enable.dsl=2

/* Modern default placement metrics. Legacy guppy metrics remain out of band. */

process GappaEDPL {
    container "${params.container__gappa}"
    label 'multithread'
    errorStrategy 'ignore'
    publishDir "${params.output}/stats", mode: 'copy'

    input:
    path dedup_jplace

    output:
    path 'edpl_list.csv', emit: edpl

    script:
    """
    set -e
    gappa examine edpl --jplace-path ${dedup_jplace} --verbose --threads ${task.cpus}
    """
}

process GappaKRD {
    container "${params.container__gappa}"
    label 'mem_veryhigh'
    errorStrategy 'ignore'
    publishDir "${params.output}/stats", mode: 'copy'

    input:
    path specimen_jplaces

    output:
    path 'krd/krd_matrix.csv.gz', emit: krd

    script:
    """
    set -e
    gappa analyze krd --jplace-path ${specimen_jplaces} --out-dir krd/ \\
      --compress --verbose --threads ${task.cpus}
    """
}

process GappaEPCA {
    container "${params.container__gappa}"
    label 'mem_veryhigh'
    errorStrategy 'ignore'
    publishDir "${params.output}/stats", mode: 'copy'

    input:
    path specimen_jplaces

    output:
    path 'epca/projection.csv', emit: projection
    path 'epca/transformation.csv', emit: transformation

    script:
    """
    set -e
    gappa analyze edgepca --jplace-path ${specimen_jplaces} --out-dir epca/ \\
      --verbose --threads ${task.cpus}
    """
}
