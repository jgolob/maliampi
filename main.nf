#!/usr/bin/env nextflow

/*
    Maliampi via nextflow!
    Steps in a typical 16S rRNA run:
    1) Make sequence variants (with dada2)
    2) create (vs load) a reference package
    3) Place SV on the reference package
    4) Classify the SV using the placements + reference package
*/

// Defaults for parameters
params.help = false
// common
params.output = '.'

// dada2-sv
params.trimLeft = 16
params.maxN = 0
params.maxEE = 'Inf'
params.truncLenF = 235
params.truncLenR = 235
params.truncQ = 2
params.errM_maxConsist = 10
params.errM_randomize = 'FALSE'
params.errM_nbases = '1e8'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run golob/maliampi <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples
                            At a minimum must have columns:
                                specimen: A unique identifier 
                                read__1: forward read
                                read__2: reverse read fq

                                optional columns:
                                batch: sequencing / library batch. Should be filename safe

                                index__1: forward index file (for checking demultiplexing)
                                index__2: reverse index file 

    Options:
      Common to all:
        --output              Directory to place outputs (default invocation dir)
                                Maliampi will create a directory structure under this directory
        -w                    Working directory. Defaults to `./work`
      SV-DADA2 options:
        --trimLeft            How far to trim on the left (default = 15)
        --maxN                (default = 0)
        --maxEE               (default = Inf)
        --truncLenF           (default = 235)
        --truncLenR           (default = 235)
        --truncQ              (default = 2)

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Load manifest!

// Check if the manifest has index files.
Channel.from(file(params.manifest))
    .splitCsv(header: true, sep: ",")
    .reduce(true){p, c ->
        return (p && (c.index__1 != null) && (c.index__2 != null))
    }.subscribe{i -> has_index = i}

// If there are index files, preceed to verifying demultiplex
if (has_index){
    Channel.from(file(params.manifest))
        .splitCsv(header: true, sep: ",")
        .map { sample -> [
            sample.specimen,
            sample.batch,
            // base names below, lacking .fq.gz
            file(sample.read__1).name.replaceAll(/.fastq.gz$/, "").replaceAll(/.fq.gz$/, ""), 
            file(sample.read__2).name.replaceAll(/.fastq.gz$/, "").replaceAll(/.fq.gz$/, ""),
            // Actual files
            file(sample.read__1),
            file(sample.read__2),
            file(sample.index__1),
            file(sample.index__2),
        ]}
        .set{ input_ch }

    // Use barcodecop to verify demultiplex
    process barcodecop {
      container "golob/barcodecop:0.4.1__bcw_0.3.0"
      label 'io_limited'
      errorStrategy "retry"

      input:
        set specimen, batch, R1_n, R2_n, file(R1), file(R2), file(I1), file(I2) from input_ch
      
      output:
        set specimen, batch, R1_n, R2_n, file("${R1_n}.bcc.fq.gz"), file("${R2_n}.bcc.fq.gz") into demultiplexed_ch
      """
      set -e

      barcodecop \
      ${I1} ${I2} \
      --match-filter \
      -f ${R1} \
      -o ${R1_n}.bcc.fq.gz &&
      barcodecop \
      ${I1} ${I2} \
      --match-filter \
      -f ${R2} \
      -o ${R2_n}.bcc.fq.gz
      """
    }
}
// Else, proceed with the file pairs as-is
else {
    Channel.from(file(params.manifest))
        .splitCsv(header: true, sep: ",")
        .map { sample -> [
            sample.specimen,
            sample.batch,
            // base names below, lacking .fq.gz
            file(sample.read__1).name.replaceAll(/.fastq.gz$/, "").replaceAll(/.fq.gz$/, ""), 
            file(sample.read__2).name.replaceAll(/.fastq.gz$/, "").replaceAll(/.fq.gz$/, ""),
            // Actual files
            file(sample.read__1),
            file(sample.read__2),
        ]}
        .set{ demultiplexed_ch }
}

// Step 1.b. next step: filter and trim reads with dada2
process dada2_ft {
    container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
    label 'io_limited'
    errorStrategy "retry"

    input:
        set specimen, batch, R1_n, R2_n, file(R1), file(R2) from demultiplexed_ch
    
    output:
        set specimen, batch, R1_n, R2_n, file("${R1_n}.dada2.ft.fq.gz"), file("${R2_n}.dada2.ft.fq.gz") into dada2_ft_ch_for_derep
        set batch, file("${R1_n}.dada2.ft.fq.gz"), file("${R2_n}.dada2.ft.fq.gz") into dada2_ft_ch_for_err
    
    """
    #!/usr/bin/env Rscript
    library('dada2'); 
    filterAndTrim(
        '${R1}', '${R1_n}.dada2.ft.fq.gz',
        '${R2}', '${R2_n}.dada2.ft.fq.gz',
        trimLeft = ${params.trimLeft},
        maxN = ${params.maxN},
        maxEE = ${params.maxEE},
        truncLen = c(${params.truncLenF}, ${params.truncLenR}),
        truncQ = ${params.truncQ},
        compress = TRUE,
        verbose = TRUE,
        multithread = ${task.cpus}
    )
    """
}

// Step 1.c. dereplicate reads with dada2

process dada2_derep {
    container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
    label 'io_limited'
    errorStrategy "retry"

    input:
        set specimen, batch, R1_n, R2_n, file(R1), file(R2) from dada2_ft_ch_for_derep
    
    output:
        set specimen, batch, R1_n, R2_n, file("${R1_n}.dada2.ft.derep.rds"), file("${R2_n}.dada2.ft.derep.rds") into dada2_derep_ch
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    derep1 <- derepFastq('${R1}');
    saveRDS(derep1, '${R1_n}.dada2.ft.derep.rds');
    derep2 <- derepFastq('${R2}');
    saveRDS(derep1, '${R2_n}.dada2.ft.derep.rds');
    """ 
}

// Step 1.d. learn errors by batch
dada2_ft_ch_for_err
    .groupTuple(by: 0)
    .set{
        dada2_ft_batches
    }

process dada2_learn_error {
    container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
    label 'multithread'
    errorStrategy "fail"
    publishDir '../working/batches'

    input:
        set val(batch), file(forwardReads), file(reverseReads) from dada2_ft_batches

    output:
        set batch, file("${batch}.R1.errM.rds"), file("${batch}.R2.errM.rds")  into dada2_batch_error_rds
        set batch, file("${batch}.R1.errM.csv"), file("${batch}.R2.errM.csv")  into dada2_batch_error_csv

    """
    #!/usr/bin/env Rscript
    library('dada2');
    errF <- learnErrors(
        strsplit(
            '${forwardReads}',
            ' '
        ),
        multithread=${task.cpus},
        MAX_CONSIST=${params.errM_maxConsist},
        randomize=${params.errM_randomize},
        nbases=${params.errM_nbases}
    );
    saveRDS(errF, "${batch}.R1.errM.rds"); 
    write.csv(errF, "${batch}.R1.errM.csv");
    errR <- learnErrors(
        strsplit(
            '${reverseReads}',
            ' '
        ),
        multithread=${task.cpus},
        MAX_CONSIST=${params.errM_maxConsist},
        randomize=${params.errM_randomize},
        nbases=${params.errM_nbases}
    );
    saveRDS(errR, "${batch}.R2.errM.rds"); 
    write.csv(errR, "${batch}.R2.errM.csv");
    """
}

