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
params.trimLeft = 0
params.maxN = 0
params.maxEE = 'Inf'
params.truncLenF = 0
params.truncLenR = 0
params.truncQ = 2
params.errM_maxConsist = 10
params.errM_randomize = 'FALSE'
params.errM_nbases = '1e8'
params.chimera_method = 'consensus'

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
        --trimLeft            How far to trim on the left (default = 0)
        --maxN                (default = 0)
        --maxEE               (default = Inf)
        --truncLenF           (default = 0)
        --truncLenR           (default = 0)
        --truncQ              (default = 2)

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.manifest == null){
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
    }.set { has_index }

// If there are index files, preceed to verifying demultiplex
if (has_index.val == true){
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
    //errorStrategy "retry"

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
        set batch, specimen, R1_n, R2_n, file("${R1_n}.dada2.ft.derep.rds"), file("${R2_n}.dada2.ft.derep.rds") into dada2_derep_ch
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    derep_1 <- derepFastq('${R1}');
    saveRDS(derep_1, '${R1_n}.dada2.ft.derep.rds');
    derep_2 <- derepFastq('${R2}');
    saveRDS(derep_2, '${R2_n}.dada2.ft.derep.rds');
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
    errorStrategy "retry"

    input:
        set val(batch), file(forwardReads), file(reverseReads) from dada2_ft_batches

    output:
        set batch, file("${batch}.R1.errM.rds"), file("${batch}.R2.errM.rds")  into dada2_batch_error_rds
        set batch, file("${batch}.R1.errM.csv"), file("${batch}.R2.errM.csv")  into dada2_batch_error_csv

    """
    #!/usr/bin/env Rscript
    library('dada2');
    errF <- learnErrors(
        unlist(strsplit(
            '${forwardReads}',
            ' ',
            fixed=TRUE
        )),
        multithread=${task.cpus},
        MAX_CONSIST=${params.errM_maxConsist},
        randomize=${params.errM_randomize},
        nbases=${params.errM_nbases}
    );
    saveRDS(errF, "${batch}.R1.errM.rds"); 
    write.csv(errF, "${batch}.R1.errM.csv");
    errR <- learnErrors(
        unlist(strsplit(
            '${reverseReads}',
            ' ',
            fixed=TRUE
        )),
        multithread=${task.cpus},
        MAX_CONSIST=${params.errM_maxConsist},
        randomize=${params.errM_randomize},
        nbases=${params.errM_nbases}
    );
    saveRDS(errR, "${batch}.R2.errM.rds"); 
    write.csv(errR, "${batch}.R2.errM.csv");
    """
}

// Step 1.e. Apply the correct error model to the dereplicated seqs
    // Join in the correct error models
    dada2_batch_error_rds
        .cross(dada2_derep_ch)
        .map{r -> [
            r[0][0], // batch
            r[1][1], // specimen
            r[1][2], // R1_n
            r[1][3], // R2_n
            file(r[1][4]), // R1
            file(r[1][5]), // R2
            file(r[0][1]), // errM_1
            file(r[0][2]), // errM_2
        ]}
        .set{ dada2_derep_w_errM_ch }


    // and use that apply the proper error model to each specimen read pair
    process dada2_dada {
        container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
        label 'multithread'
        errorStrategy "retry"

        publishDir '../working/dada2/dada/', mode: 'copy'

        input:
            set batch, specimen, R1_n, R2_n, file(R1), file(R2), file(errM_1), file(errM_2) from dada2_derep_w_errM_ch

        output:
            set batch, specimen, file(R1), file(R2), file("${R1_n}.dada2.dada.rds"), file("${R2_n}.dada2.dada.rds") into dada2_dada_sp_ch

        """
        #!/usr/bin/env Rscript
        library('dada2');
        errM_1 <- readRDS('${errM_1}');
        derep_1 <- readRDS('${R1}');
        dadaResult_1 <- dada(derep_1, err=errM_1, multithread=${task.cpus}, verbose=FALSE);
        saveRDS(dadaResult_1, '${R1_n}.dada2.dada.rds');
        errM_2 <- readRDS('${errM_2}');
        derep_2 <- readRDS('${R2}');
        dadaResult_2 <- dada(derep_2, err=errM_2, multithread=${task.cpus}, verbose=FALSE);
        saveRDS(dadaResult_2, '${R2_n}.dada2.dada.rds');
        """
    }

// Step 1.f. Merge reads, using the dereplicated seqs and applied model

process dada2_merge {
        container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
        label 'multithread'
        //errorStrategy "retry"

        publishDir '../working/dada2/merged/', mode: 'copy'

        input:
            set batch, specimen, file(R1), file(R2), file(R1dada), file(R2dada) from dada2_dada_sp_ch

        output:
            set batch, specimen, file("${specimen}.dada2.merged.rds") into dada2_sp_merge

        """
        #!/usr/bin/env Rscript
        library('dada2');
        dada_1 <- readRDS('${R1dada}');
        derep_1 <- readRDS('${R1}');
        dada_2 <- readRDS('${R2dada}');
        derep_2 <- readRDS('${R2}');        
        merger <- mergePairs(
            dada_1, derep_1,
            dada_2, derep_2,
            verbose=TRUE,
        );
        saveRDS(merger, "${specimen}.dada2.merged.rds");
        """
    }

// Step 1.g. Make a seqtab

process dada2_seqtab_sp {
        container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
        label 'io_limited'
        errorStrategy "retry"

        input:
            set batch, specimen, file(merged) from dada2_sp_merge

        output:
            set batch, specimen, file("${specimen}.dada2.seqtab.rds") into dada2_sp_seqtab

        """
        #!/usr/bin/env Rscript
        library('dada2');
        merged <- readRDS('${merged}');
        seqtab <- makeSequenceTable(merged);
        rownames(seqtab) <- c('${specimen}');
        saveRDS(seqtab, '${specimen}.dada2.seqtab.rds');
        """
    }

// Step 1.h. Combine seqtabs
    // Do this by batch to help with massive data sets
    dada2_sp_seqtab
        .groupTuple(by: 0)
        .set{ dada2_batch_seqtabs_ch }

    process dada2_seqtab_batch_combine {
        container 'golob/dada2-fast-combineseqtab:0.2.0_BCW_0.30A'
        label 'io_limited'
        //errorStrategy "retry"

        input:
            set val(batch), val(specimens), file(sp_seqtabs_rds) from dada2_batch_seqtabs_ch

        output:
            set batch, file("${batch}.dada2.seqtabs.rds") into dada2_batch_seqtab_ch

        """
        combine_seqtab \
        --rds ${batch}.dada2.seqtabs.rds \
        --seqtabs ${sp_seqtabs_rds}
        """
    }
    // Then combine the batch seqtabs into the one seqtab to rule them all
    dada2_batch_seqtab_ch
        .collect{ it[1] }
        .set{ batch_seqtab_files }



    process dada2_seqtab_combine_all {
        container 'golob/dada2-fast-combineseqtab:0.2.0_BCW_0.30A'
        label 'io_limited'
        // errorStrategy "retry"

        input:
            file(batch_seqtab_files)

        output:
            file("combined.dada2.seqtabs.rds") into combined_seqtab

        """
        combine_seqtab \
        --rds combined.dada2.seqtabs.rds \
        --seqtabs ${batch_seqtab_files}
        """
    }
// Step 1.i. Remove chimera on combined seqtab
    process dada2_remove_bimera {
        container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
        label 'multithread'
        publishDir "${params.output}/sv/", mode: 'copy'

        input:
            file(combined_seqtab)

        output:
            file("dada2.combined.seqtabs.nochimera.rds") into final_seqtab_rds
            file("dada2.combined.seqtabs.nochimera.csv") into final_seqtab_csv

        """
        #!/usr/bin/env Rscript
        library('dada2');
        seqtab <- readRDS('${combined_seqtab}');
        seqtab_nochim <- removeBimeraDenovo(
            seqtab,
            method = '${params.chimera_method}',
            multithread = ${task.cpus}
        );
        saveRDS(seqtab_nochim, 'dada2.combined.seqtabs.nochimera.rds'); 
        write.csv(seqtab_nochim, 'dada2.combined.seqtabs.nochimera.csv', na='');
        print((sum(seqtab) - sum(seqtab_nochim)) / sum(seqtab));
        """
    }
// Step 1.j. Transform output to be pplacer and mothur style
    process dada2_convert_output {
        container 'golob/dada2-pplacer:0.4.1__bcw_0.3.1'
        label 'io_limited'
        publishDir "${params.output}/sv/", mode: 'copy'

        input:
            file(final_seqtab_csv)

        output:
            file("dada2.sv.fasta") into dada2_sv_fasta
            file("dada2.sv.map.csv") into dada2_sv_map
            file("dada2.sv.weights.csv") into dada2_sv_weights
            file("dada2.sv.shared.txt") into dada2_sv_sharetable


        """
        dada2-seqtab-to-pplacer \
        -s ${final_seqtab_csv} \
        -f dada2.sv.fasta \
        -m dada2.sv.map.csv \
        -w dada2.sv.weights.csv \
        -t dada2.sv.shared.txt
        """
    }

// */

