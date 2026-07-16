nextflow.enable.dsl=2

process GoodsFilter {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'mem_medium'
    publishDir "${params.output}/goods", mode: 'copy'

    input:
    path sv_h5ad

    output:
    path 'sv.goods.h5ad', emit: sv_h5ad
    path 'goods.convergence.parquet', emit: convergence
    path 'goods.curves.parquet', emit: curves

    script:
    keep_nonconverged = params.goods_keep_nonconverged ? '--keep-nonconverged' : ''
    """
    set -euo pipefail
    maliampi-goods-filter ${sv_h5ad} sv.goods.h5ad \
        --convergence-parquet goods.convergence.parquet \
        --curves-parquet goods.curves.parquet \
        --convergence-delta ${params.goods_convergence_delta} \
        --min-reads ${params.goods_min_reads} \
        --min-prevalence ${params.goods_min_prevalence} \
        --seed ${params.goods_seed} ${keep_nonconverged}
    """
}
