container__barcodecop = "golob/barcodecop:0.4.1__bcw_0.3.0"

workflow preprocess_wf {
    take: indexed_ch
    take: paired_ch

    main:

    // Indexed rows into a channel for barcodecop
    indexed_ch.map{ sample -> [
        sample.specimen,
        sample.batch,
        file(sample.R1),
        file(sample.R2),
        file(sample.I1),
        file(sample.I2),
    ]}.set{ to_bcc_ch }

    // Run barcodecop
    barcodecop(to_bcc_ch)

    // Raise an error if any of the samples fail barcodecop
    bcc_results = barcodecop.out.branch {
        empty: (file(it[2]).isEmpty() || file(it[3]).isEmpty())
        valid: true
    }

    // Mix together the reads with no index with the reads with verified demultiplex
     paired_ch
        .map { sample -> [
            sample.specimen,
            sample.batch,
            file(sample.R1),
            file(sample.R2),
        ]}
        .mix(bcc_results.valid)
        .set{ demultiplexed_ch }
    emit:
        valid = demultiplexed_ch
        empty = bcc_results.empty
}

// Use barcodecop to verify demultiplex
process barcodecop {
    container "${container__barcodecop}"
    label 'io_limited'
    errorStrategy "retry"

    input:
    tuple val(specimen), val(batch), file(R1), file(R2), file(I1), file(I2)
    
    output:
    tuple val(specimen), val(batch), file("${R1.getSimpleName()}.bcc.fq.gz"), file("${R2.getSimpleName()}.bcc.fq.gz")

    """
    set -e

    barcodecop \
    ${I1} ${I2} \
    --match-filter \
    -f ${R1} \
    -o ${R1.getSimpleName()}.bcc.fq.gz &&
    barcodecop \
    ${I1} ${I2} \
    --match-filter \
    -f ${R2} \
    -o ${R2.getSimpleName()}.bcc.fq.gz
    """
}


