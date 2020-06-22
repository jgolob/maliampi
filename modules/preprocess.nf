container__barcodecop = "golob/barcodecop:0.4.1__bcw_0.3.0"

workflow preprocess_wf {
    take: indexed_ch
    take: paired_ch
    take: unpaired_ch

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

    // For the unpaired reads, figure out if they are miseq or pyro using the 'data_type' column. 
    // If non-existant, assume it is miseq
    unpaired_ch
        .branch {
            pyro: (it.data_type == '16S_pyro') | (it.data_type == 'pyro')
            miseq: true
        }.set { single_end_branch_ch }

    single_end_branch_ch.miseq.map{ sample -> [
        sample.specimen,
        sample.batch,
        file(sample.R1)
    ]}.set {
        miseq_se
    }

    single_end_branch_ch.pyro.map{ sample -> [
        sample.specimen,
        sample.batch,
        file(sample.R1)
    ]}.set {
        pyro_ch
    }    

    // Final outputs
    emit:
        miseq_pe = demultiplexed_ch
        miseq_se = miseq_se
        pyro = pyro_ch
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

process output_failed {
    container "${container__barcodecop}"
    label 'io_limited'
    publishDir "${params.output}/sv/", mode: 'copy'
    errorStrategy 'retry'

    input:
        tuple val(specimens), val(reasons)
    output:
        file ("failed_specimens.csv")

    """
    #!/usr/bin/env python
    import csv
    import re
    specimens = re.sub(r'\\[|\\]', "", "${specimens}").split(',')
    reasons = re.sub(r'\\[|\\]', "", "${reasons}").split(',')

    with open("failed_specimens.csv", 'wt') as out_h:
        w = csv.writer(out_h)
        w.writerow([
            'specimen',
            'reason'
        ])
        for sp, reason in zip(specimens, reasons):
            w.writerow([sp.strip(), reason.strip()])
    """
}