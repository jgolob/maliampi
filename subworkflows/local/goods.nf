#!/usr/bin/env nextflow

/* Opt-in standalone Good's-coverage filtering for canonical SV H5AD. */
nextflow.enable.dsl=2

include { GoodsFilter } from '../../modules/goods'

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
    if (params.sv_h5ad == null) {
        error "Good's filtering requires --sv_h5ad"
    }
    goods(channel.fromPath(params.sv_h5ad, checkIfExists: true))
}
