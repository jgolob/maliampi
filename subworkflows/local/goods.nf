#!/usr/bin/env nextflow

/* Opt-in standalone Good's-coverage filtering for canonical SV H5AD. */

params.help = false
params.sv_h5ad = null
params.output = '.'

nextflow.enable.dsl=2

include { GoodsFilter } from '../../modules/goods'

def helpMessage() {
    log.info """
    ─────────────────────────────────────
    maliampi / goods
    ─────────────────────────────────────
    Good's-coverage filtering for canonical SV H5AD artifacts.
    Removes specimens that have not reached sequencing saturation.

    Usage:
      nextflow run jgolob/maliampi/subworkflows/local/goods.nf [options]

    Required:
      --sv_h5ad       H5AD file with sequence variants

    Options:
      --output        Output directory (default: .)
      --help          Show this help message
    """.stripIndent()
}

workflow goods {
    take:
    sv_h5ad

    main:
    GoodsFilter(sv_h5ad)

    emit:
    sv_h5ad = GoodsFilter.out.sv_h5ad
    convergence = GoodsFilter.out.convergence
    curves = GoodsFilter.out.curves
}

workflow {
    if (params.help || params.sv_h5ad == null) {
        helpMessage()
        if (!params.help) {
            error "Missing required parameter(s): --sv_h5ad"
        }
        return
    }
    goods(channel.fromPath(params.sv_h5ad, checkIfExists: true))
}
