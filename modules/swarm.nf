#!/usr/bin/env nextflow

/* Standalone Swarm SV workflow. Its internal and published contract is H5AD. */
nextflow.enable.dsl=2

include { read_manifest } from './manifest'
include { output_failed; preprocess_wf } from './preprocess'

workflow swarm_wf {
    take:
    miseq_pe_ch
    miseq_se_ch
    pyro_ch

    main:
    MergePairs(miseq_pe_ch)

    reads = MergePairs.out
        .mix(miseq_se_ch)
        .mix(pyro_ch)
    FilterAndTrim(reads)
    SpecimenDereplicate(FilterAndTrim.out)
    PrepareSpecimenRegistry(SpecimenDereplicate.out)

    registries = PrepareSpecimenRegistry.out.registry
        .map { specimen, batch, registry -> registry }
        .collect()
    CombineDsvRegistries(registries)
    SwarmCluster(CombineDsvRegistries.out.dsv_fasta)
    FilterSwarmSeeds(SwarmCluster.out.seeds)
    ChimeraRemoval(FilterSwarmSeeds.out)
    FinalizeSwarm(
        CombineDsvRegistries.out.registry,
        SwarmCluster.out.clusters,
        ChimeraRemoval.out,
    )

    emit:
    sv_h5ad = FinalizeSwarm.out.sv_h5ad
    stats = FinalizeSwarm.out.stats
}

process MergePairs {
    container "${params.container__vsearch}"
    label 'io_limited'
    errorStrategy 'ignore'

    input:
    tuple val(specimen), val(batch), path(R1), path(R2)

    output:
    tuple val(specimen), val(batch), path('merged.fastq.gz')

    script:
    """
    set -euo pipefail
    vsearch --fastq_mergepairs ${R1} --reverse ${R2} \
        --fastqout merged.fastq --threads ${task.cpus} --fastq_eeout
    gzip merged.fastq
    """
}

process FilterAndTrim {
    container "${params.container__vsearch}"
    label 'io_limited'
    errorStrategy 'ignore'

    input:
    tuple val(specimen), val(batch), path(reads)

    output:
    tuple val(specimen), val(batch), path('filtered.fasta.gz')

    script:
    """
    set -euo pipefail
    vsearch --fastq_filter ${reads} \
        --fastq_maxee ${params.maxEE} \
        --fastq_maxns ${params.maxN} \
        --fastq_truncqual ${params.truncQ} \
        --fastq_stripleft ${params.trimLeft} \
        --fastq_stripright ${params.truncLenR} \
        --threads ${task.cpus} \
        --fastaout filtered.fasta
    gzip filtered.fasta
    """
}

process SpecimenDereplicate {
    container "${params.container__vsearch}"
    label 'multithread'
    errorStrategy 'ignore'

    input:
    tuple val(specimen), val(batch), path(reads)

    output:
    tuple val(specimen), val(batch), path('derep.fasta.gz')

    script:
    """
    set -euo pipefail
    vsearch --derep_fulllength ${reads} \
        --strand plus --sizeout --fasta_width 0 --threads ${task.cpus} \
        --output derep.fasta
    gzip derep.fasta
    """
}

process PrepareSpecimenRegistry {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'
    errorStrategy 'ignore'

    input:
    tuple val(specimen), val(batch), path(derep_fasta)

    output:
    tuple val(specimen), val(batch), path('specimen.*.registry.parquet'), emit: registry

    script:
    specimen_arg = specimen.toString().replace("'", "'\"'\"'")
    batch_arg = (batch ?: '').toString().replace("'", "'\"'\"'")
    """
    set -euo pipefail
    maliampi-swarm-specimen ${derep_fasta} specimen.${task.index}.registry.parquet \
        --specimen '${specimen_arg}' --batch '${batch_arg}'
    """
}

process CombineDsvRegistries {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'mem_medium'

    input:
    path registries

    output:
    path 'dsv.fasta.gz', emit: dsv_fasta
    path 'dsv.registry.parquet', emit: registry

    script:
    """
    set -euo pipefail
    maliampi-swarm-combine --output-fasta dsv.fasta.gz \
        --output-registry dsv.registry.parquet ${registries}
    """
}

process SwarmCluster {
    container "${params.container__swarm}"
    label 'multithread'

    input:
    path dsv_fasta

    output:
    path 'swarm.clusters', emit: clusters
    path 'swarm.seeds.fasta.gz', emit: seeds

    script:
    """
    set -euo pipefail
    gzip -dc ${dsv_fasta} > dsv.fasta
    swarm -d ${params.swarm_d} -f -z -t ${task.cpus} \
        -o swarm.clusters -w swarm.seeds.fasta dsv.fasta
    gzip swarm.seeds.fasta
    """
}

process FilterSwarmSeeds {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path seeds

    output:
    path 'swarm.seeds.nonsingleton.fasta.gz'

    script:
    """
    set -euo pipefail
    maliampi-swarm-filter-seeds ${seeds} swarm.seeds.nonsingleton.fasta.gz
    """
}

process ChimeraRemoval {
    container "${params.container__vsearch}"
    label 'multithread'

    input:
    path seeds

    output:
    path 'swarm.nonchimera.fasta.gz'

    script:
    """
    set -euo pipefail
    gzip -dc ${seeds} > seeds.fasta
    vsearch --uchime_denovo seeds.fasta \
        --nonchimeras swarm.nonchimera.fasta --threads ${task.cpus}
    gzip swarm.nonchimera.fasta
    """
}

process FinalizeSwarm {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'mem_medium'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
    path registry
    path clusters
    path nonchimera

    output:
    path 'sv.h5ad', emit: sv_h5ad
    path 'swarm.stats.parquet', emit: stats

    script:
    project_arg = params.project_id.toString().replace("'", "'\"'\"'")
    dataset_arg = params.dataset_id.toString().replace("'", "'\"'\"'")
    """
    set -euo pipefail
    maliampi-swarm-finalize ${registry} ${clusters} ${nonchimera} \
        sv.h5ad swarm.stats.parquet \
        --project-id '${project_arg}' --dataset-id '${dataset_arg}'
    """
}

workflow {
    if (params.manifest == null) {
        error 'Swarm requires --manifest'
    }
    if (params.project_id == null) {
        error 'Swarm requires --project_id'
    }
    if (params.dataset_id == null) {
        error 'Swarm requires --dataset_id'
    }

    manifest = read_manifest(channel.fromPath(params.manifest, checkIfExists: true))
    preprocess_wf(
        manifest.valid_paired_indexed,
        manifest.valid_paired,
        manifest.valid_unpaired,
    )
    swarm_wf(
        preprocess_wf.out.miseq_pe,
        preprocess_wf.out.miseq_se,
        preprocess_wf.out.pyro,
    )

    failures = manifest.other.map { row -> [row.specimen, 'failed at manifest'] }
        .mix(preprocess_wf.out.empty.map { row -> [row[0], 'preprocessing'] })
        .collect()
        .map { rows -> [
            rows.collect { row -> row[0] },
            rows.collect { row -> row[1] },
        ] }
    output_failed(failures)
}
