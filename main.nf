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
        --taxdmp              Path to taxdmp.zip. If not provided, it will be downloaded
      Ref Package options:
        --repo_fasta          FASTA file containing reference sequences

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.manifest == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

//
//  START STEP 1: Sequence variants

// Create a channel for specimens whose reads are filtered away at one of the steps
// The structure is a tuple, first is the specimen ID, the second is the step.
filtered_specimens_ch = Channel.create()

// Load manifest!

// For each, figure out if an index is available, and split into with index and without channels
// Also check that we have at least a specimen and read__1 and read__2 provided

input_w_index_ch = Channel.create()
input_no_index_ch = Channel.create()
input_invalid_ch = Channel.create()

Channel.from(file(params.manifest))
    .splitCsv(header: true, sep: ",")
    .choice(
        input_invalid_ch,
        input_no_index_ch,
        input_w_index_ch,
    ) {
        r -> if (
            (r.specimen == null) ||
            (r.read__1 == null) ||
            (r.read__2 == null) ||
            (r.specimen == "") ||
            (r.read__1 == "") ||
            (r.read__2 == "")
        ) return 0;
        else if (
            (r.index__1 != null) && 
            (r.index__2 != null) && 
            (r.index__1 != "") &&
            (r.index__2 != "")
            ) return 2;
        else return 1;
    }

input_invalid_ch.subscribe{
    print "Missing required specimen, read__1 or read__2: "
    println it
}


// Now check to be sure the files exist and are not empty.
// If any file in a row is empty, make a note of it and proceed with the remainder

input_w_index_invalid_ch = Channel.create()
input_w_index_valid_ch = Channel.create()
input_w_index_ch.choice(
    input_w_index_invalid_ch,
    input_w_index_valid_ch
) {
    r -> (file(r.read__1).isEmpty() || file(r.read__2).isEmpty() || file(r.index__1).isEmpty() || file(r.index__2).isEmpty()) ? 0 : 1
}

input_no_index_invalid_ch = Channel.create()
input_no_index_valid_ch = Channel.create()
input_no_index_ch.choice(
    input_no_index_invalid_ch,
    input_no_index_valid_ch
) {
    r -> (file(r.read__1).isEmpty() || file(r.read__2).isEmpty()) ? 0 : 1
}

// Handle the invalids. Here we are just going to print them as we see them.
input_no_index_invalid_ch.subscribe{
    print "Missing / Empty files for: "
    println it
}
input_w_index_invalid_ch.subscribe{
    print "Missing / Empty files for: "
    println it
}

// For those with an index, make a channel for barcodecop
input_w_index_valid_ch
    .map{ sample -> [
        sample.specimen,
        sample.batch,
        file(sample.read__1),
        file(sample.read__2),
        file(sample.index__1),
        file(sample.index__2),
    ]}
    .set{ to_bcc_ch }

// Use barcodecop to verify demultiplex
process barcodecop {
    container "golob/barcodecop:0.4.1__bcw_0.3.0"
    label 'io_limited'
    errorStrategy "retry"

    input:
    set specimen, batch, file(R1), file(R2), file(I1), file(I2) from to_bcc_ch
    
    output:
    set specimen, batch, file("${R1.getSimpleName()}.bcc.fq.gz"), file("${R2.getSimpleName()}.bcc.fq.gz") into from_bcc_ch
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


bcc_to_ft_ch = Channel.create()
bcc_empty_ch = Channel.create()
from_bcc_ch
    .choice(
        bcc_to_ft_ch,
        bcc_empty_ch
    ) {
        r -> (file(r[2]).isEmpty() || file(r[3]).isEmpty()) ? 1 : 0
    }

bcc_empty_ch.subscribe {
    filtered_specimens_ch.bind([it[0], 'BCC']);
    print "No reads from "
    print it[0]
    println " survived BCC"
}

// Else, proceed with the file pairs as-is

input_no_index_valid_ch
    .map { sample -> [
        sample.specimen,
        sample.batch,
        // Actual files
        file(sample.read__1),
        file(sample.read__2),
    ]}
    .set{ no_index_to_ft_ch }

no_index_to_ft_ch.mix(bcc_to_ft_ch).set{
    demultiplexed_ch
}

// Step 1.b. next step: filter and trim reads with dada2
process dada2_ft {
    container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
    label 'io_limited'
    errorStrategy "retry"

    input:
        set specimen, batch, file(R1), file(R2) from demultiplexed_ch
    
    output:
        set specimen, batch, file("${R1.getSimpleName()}.dada2.ft.fq.gz"), file("${R2.getSimpleName()}.dada2.ft.fq.gz") into dada2_post_ft_ch
    
    """
    #!/usr/bin/env Rscript
    library('dada2'); 
    filterAndTrim(
        '${R1}', '${R1.getSimpleName()}.dada2.ft.fq.gz',
        '${R2}', '${R2.getSimpleName()}.dada2.ft.fq.gz',
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

// Filter out specimens with no surviving reads after FT
dada2_ft_ch = Channel.create()
dada2_ft_empty_ch = Channel.create()

dada2_post_ft_ch
    .choice(
        dada2_ft_ch,
        dada2_ft_empty_ch
    ) {
        ( file(it[2]).isEmpty() || file(it[3]).isEmpty() ) ? 1 : 0
    }

dada2_ft_ch.into{
    dada2_ft_for_derep_ch;
    dada2_ft_ch_for_err
}

dada2_ft_empty_ch.subscribe{
    filtered_specimens_ch.bind([it[0], 'FT']);
    print "No reads from "
    print it[0]
    println " survived filter/trim"    
}

// Step 1.c. dereplicate reads with dada2

process dada2_derep {
    container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
    label 'io_limited'
    errorStrategy "retry"

    input:
        set specimen, batch, file(R1), file(R2) from dada2_ft_for_derep_ch
    
    output:
        set batch, specimen, file("${R1.getSimpleName()}.dada2.ft.derep.rds"), file("${R2.getSimpleName()}.dada2.ft.derep.rds") into dada2_derep_ch
    
    """
    #!/usr/bin/env Rscript
    library('dada2');
    derep_1 <- derepFastq('${R1}');
    saveRDS(derep_1, '${R1.getSimpleName()}.dada2.ft.derep.rds');
    derep_2 <- derepFastq('${R2}');
    saveRDS(derep_2, '${R2.getSimpleName()}.dada2.ft.derep.rds');
    """ 
}

// Step 1.d. learn errors by batch
dada2_ft_ch_for_err
    .groupTuple(by: 1)
    .set{
        dada2_ft_batches
    }

process dada2_learn_error {
    container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
    label 'multithread'
    errorStrategy "retry"
    //maxForks 5

    input:
        set val(batch_specimens), val(batch), file(forwardReads), file(reverseReads) from dada2_ft_batches

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
            file(r[1][2]), // R1
            file(r[1][3]), // R2
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
            set batch, specimen, file(R1), file(R2), file(errM_1), file(errM_2) from dada2_derep_w_errM_ch

        output:
            set batch, specimen, file(R1), file(R2), file("${R1.getSimpleName()}.dada2.dada.rds"), file("${R2.getSimpleName()}.dada2.dada.rds") into dada2_dada_sp_ch

        """
        #!/usr/bin/env Rscript
        library('dada2');
        errM_1 <- readRDS('${errM_1}');
        derep_1 <- readRDS('${R1}');
        dadaResult_1 <- dada(derep_1, err=errM_1, multithread=${task.cpus}, verbose=FALSE);
        saveRDS(dadaResult_1, '${R1.getSimpleName()}.dada2.dada.rds');
        errM_2 <- readRDS('${errM_2}');
        derep_2 <- readRDS('${R2}');
        dadaResult_2 <- dada(derep_2, err=errM_2, multithread=${task.cpus}, verbose=FALSE);
        saveRDS(dadaResult_2, '${R2.getSimpleName()}.dada2.dada.rds');
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
            set batch, specimen, file("${specimen}.dada2.merged.rds") into dada2_sp_post_merge_ch

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

dada2_sp_merge_ch = Channel.create()
dada2_sp_merge_empty_ch = Channel.create()

dada2_sp_post_merge_ch
    .choice(
        dada2_sp_merge_ch,
        dada2_sp_merge_empty_ch 
    ) {
        file(it[2]).isEmpty() ? 1 : 0
    }

dada2_sp_merge_empty_ch
    .subscribe{
        filtered_specimens_ch.bind([it[1], 'merge']);
        print "No reads from "
        print it[1]
        println " survived merge"         
    }

// Step 1.g. Make a seqtab

process dada2_seqtab_sp {
        container 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'
        label 'io_limited'
        errorStrategy "retry"

        input:
            set batch, specimen, file(merged) from dada2_sp_merge_ch

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
            file "combined.dada2.seqtabs.rds" into combined_seqtab

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
            file "dada2.combined.seqtabs.nochimera.rds" into final_seqtab_rds
            file "dada2.combined.seqtabs.nochimera.csv" into final_seqtab_csv

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
            file "dada2.sv.fasta" into dada2_sv_fasta_f
            file "dada2.sv.fasta"  into sv_fasta_f
            file "dada2.sv.map.csv"  into sv_map_f
            file "dada2.sv.weights.csv" into sv_weights_f
            file "dada2.sv.shared.txt" into dada2_sv_sharetable_f


        """
        dada2-seqtab-to-pplacer \
        -s ${final_seqtab_csv} \
        -f dada2.sv.fasta \
        -m dada2.sv.map.csv \
        -w dada2.sv.weights.csv \
        -t dada2.sv.shared.txt
        """
    }

// Step 1.k. Output any failed specimens, and the step at which they failed.
    filtered_specimens_ch.reduce('specimen,step\n'){
        p, c -> return p+="${c[0]},${c[1]}\n";
    }.set {filtered_specimens_str }
    process output_failed {
        container = "golob/dada2-pplacer:0.4.1__bcw_0.3.1"
        publishDir "${params.output}/sv/", mode: 'copy'

        input:
            val filtered_specimens_str
        output:
            file "failed_specimens.csv"

        """
        echo "${filtered_specimens_str}" > failed_specimens.csv
        """

    }

//
//  END STEP 1: Sequence variants 
//


//
//  START STEP 2: Reference package
//
params.repo_min_id = 0.8
params.repo_max_accepts = 10
params.cmalign_mxsize = 8196
params.raxml_model = 'GTRGAMMA'
params.raxml_parsiomony_seed = 12345


// Step 2.a. Search the repo for matches
// load the repo
Channel.value(file(params.repo_fasta))
    .set{ repo_fasta}
Channel.value(file(params.repo_si))
    .set{ repo_si }

process refpkgSearchRepo {
    container 'golob/vsearch:2.7.1_bcw_0.2.0'
    label = 'multithread'

    input:
        file sv_fasta_f
        file repo_fasta
    
    output:
        file "${repo_fasta}.uc" into sv_repo_uc_f
        file "${repo_fasta}.sv.nohit.fasta" into sv_repo_nohit_f
        file "${repo_fasta}.repo.recruits.fasta" into repo_recruits_f

    """
    vsearch \
    --threads=${task.cpus} \
    --usearch_global ${sv_fasta_f} \
    --db ${repo_fasta} \
    --id=${params.repo_min_id} \
    --strand both \
    --uc=${repo_fasta}.uc --uc_allhits \
    --notmatched=${repo_fasta}.sv.nohit.fasta \
    --dbmatched=${repo_fasta}.repo.recruits.fasta \
    --maxaccepts=${params.repo_max_accepts} \
    | tee -a vsearch.log
    """
}


// Step 2.xx Filter SeqInfo to recruits
process filterSeqInfo {
    container = 'golob/fastatools:0.7.1__bcw.0.3.1'
    label = 'io_limited'

    input:
        file repo_recruits_f
        file repo_si
    
    output:
        file 'refpkg.seq_info.csv' into refpkg_si_f

    """
    #!/usr/bin/env python
    import fastalite
    import csv

    with open('${repo_recruits_f}', 'rt') as fasta_in:
        seq_ids = {sr.id for sr in fastalite.fastalite(fasta_in)}
    with open('${repo_si}', 'rt') as si_in, open('refpkg.seq_info.csv', 'wt') as si_out:
        si_reader = csv.DictReader(si_in)
        si_writer = csv.DictWriter(si_out, si_reader.fieldnames)
        si_writer.writeheader()
        for r in si_reader:
            if r['seqname'] in seq_ids:
                si_writer.writerow(r)
    """
}

// Step 2.xx (get) or build a taxonomy db

if ( (params.taxdmp == null) || file(params.taxdmp).isEmpty() ) {
    process DlBuildTaxtasticDB {
        container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
        label = 'io_limited'
        // errorStrategy = 'retry'

        output:
            file "taxonomy.db" into taxonomy_db_f

        afterScript "rm -rf dl/"


        """
        mkdir -p dl/ && \
        taxit new_database taxonomy.db -p dl/
        """

    }

} else {
    taxdump_zip_f = file(params.taxdmp)
    process buildTaxtasticDB {
        container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
        label = 'io_limited'
        // errorStrategy = 'retry'

        input:
            file taxdump_zip_f

        output:
            file "taxonomy.db" into taxonomy_db_f

        """
        taxit new_database taxonomy.db -z ${taxdump_zip_f}
        """
    }
}


// Step 2.xx Confirm seq info taxonomy matches taxdb

process confirmSI {
    container = 'golob/seqinfo_taxonomy_sync:0.2.1__bcw.0.3.0'
    label = 'io_limited'

    input:
        file taxonomy_db_f
        file refpkg_si_f
    
    output:
        file "${refpkg_si_f.baseName}.corr.csv" into refpkg_si_corr_f
    
    """
    seqinfo_taxonomy_sync.py \
    ${refpkg_si_f} ${refpkg_si_f.baseName}.corr.csv \
    --db ${taxonomy_db_f} --email ${params.email}
    """


}

// Step 2.xx Align recruited seqs

process alignRepoRecruits {
    container = 'golob/infernal:1.1.2_bcw_0.2.0'
    label = 'mem_veryhigh'

    input:
        file repo_recruits_f
    
    output:
        file "recruits.aln.scores" into recruit_aln_scores_f
        file "recruits.aln.sto" into recruits_aln_sto_f
    
    """
    cmalign \
    --cpu ${task.cpus} --noprob --dnaout --mxsize ${params.cmalign_mxsize} \
    --sfile recruits.aln.scores -o recruits.aln.sto \
    /cmalign/data/SSU_rRNA_bacteria.cm ${repo_recruits_f}
    """
}

// Step 2.xx Convert alignment from STO -> FASTA format
process convertAlnToFasta {
    container = 'golob/fastatools:0.7.1__bcw.0.3.1'
    label = 'io_limited'
    errorStrategy "retry"

    input: 
        file recruits_aln_sto_f
    
    output:
        file "recruits.aln.fasta" into recruits_aln_fasta_f
    
    """
    #!/usr/bin/env python
    from Bio import AlignIO

    with open('recruits.aln.fasta', 'wt') as out_h:
        AlignIO.write(
            AlignIO.read(
                open('${recruits_aln_sto_f}', 'rt'),
                'stockholm'
            ),
            out_h,
            'fasta'
        )
    """
}
// Step 2.xx Make a tree from the alignment.
process raxmlTree {
    container = 'golob/raxml:8.2.11_bcw_0.3.0'
    label = 'mem_veryhigh'
    errorStrategy = 'retry'

    input:
        file recruits_aln_fasta_f
    
    output:
        file "RAxML_bestTree.refpkg" into refpkg_tree_f
        file "RAxML_info.refpkg" into refpkg_tree_stats_f
    
    """
    raxml \
    -n refpkg \
    -m ${params.raxml_model} \
    -s ${recruits_aln_fasta_f} \
    -p ${params.raxml_parsiomony_seed} \
    -T ${task.cpus} && \
    ls -l -h 
    """
}



// Step 2.xx Make a tax table for the refpkg sequences
process taxtableForSI {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'
    // errorStrategy = 'retry'

    input:
        file taxonomy_db_f 
        file refpkg_si_f 
    output:
        file "refpkg.taxtable.csv" into refpkg_tt_f

    """
    taxit taxtable ${taxonomy_db_f} \
    --seq-info ${refpkg_si_f} \
    --outfile refpkg.taxtable.csv
    """
}

// Step 2.xx Obtain the CM used for the alignment
process obtainCM {
    container = 'golob/infernal:1.1.2_bcw_0.2.0'
    label = 'io_limited'

    output:
        file "refpkg.cm" into refpkg_cm
    
    """
    cp /cmalign/data/SSU_rRNA_bacteria.cm refpkg.cm
    """
}

// Step 2.xx Combine into a refpkg
process combineRefpkg {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'

    afterScript("rm -rf refpkg/*")
    publishDir "${params.output}/refpkg/", mode: 'copy'

    input:
        file recruits_aln_fasta_f
        file recruits_aln_sto_f
        file refpkg_tree_f 
        file refpkg_tree_stats_f 
        file refpkg_tt_f
        file refpkg_si_corr_f
        file refpkg_cm
    
    output:
        file "refpkg.tgz" into refpkg_tgz_f
    
    """
    taxit create --locus 16S \
    --package-name refpkg \
    --clobber \
    --aln-fasta ${recruits_aln_fasta_f} \
    --aln-sto ${recruits_aln_sto_f} \
    --tree-file ${refpkg_tree_f} \
    --tree-stats ${refpkg_tree_stats_f} \
    --taxonomy ${refpkg_tt_f} \
    --seq-info ${refpkg_si_corr_f} \
    --profile ${refpkg_cm} && \
    ls -l refpkg/ && \
    tar czvf refpkg.tgz  -C refpkg/ .
    """
}

//
//  END STEP 2: Reference package
//

//
//  START STEP 3: Placement
//
params.pplacer_prior_lower = 0.01

//  Step 3.a. Align SV

process alignSV {
    container = 'golob/infernal:1.1.2_bcw_0.2.0'
    label = 'mem_veryhigh'

    input:
        file sv_fasta_f
    
    output:
        file "sv.aln.scores" into sv_aln_scores_f
        file "sv.aln.sto" into sv_aln_sto_f
    
    """
    cmalign \
    --cpu ${task.cpus} --noprob --dnaout --mxsize ${params.cmalign_mxsize} \
    --sfile sv.aln.scores -o sv.aln.sto \
    /cmalign/data/SSU_rRNA_bacteria.cm ${sv_fasta_f}
    """
}

//  Step 3.b. Combine SV and refpkg alignment
// First extract the alignment from the refpkg

process extractRefpkgAln {
    container = "golob/fastatools:0.7.1__bcw.0.3.1"
    label = 'io_limited'

    input:
        file refpkg_tgz_f
    
    output:
        file "refpkg.aln.fasta"
        file "refpkg.aln.sto" into refpkg_aln_sto_f
        
    """
    #!/usr/bin/env python

    import tarfile
    import json
    from Bio import AlignIO
    import os

    tar_h = tarfile.open('${refpkg_tgz_f}')
    tar_contents_dict = {os.path.basename(f.name): f for f in tar_h.getmembers()}
    print(tar_contents_dict)
    contents = json.loads(
        tar_h.extractfile(
            tar_contents_dict['CONTENTS.json']
        ).read().decode('utf-8')
    )
    aln_fasta_intgz = contents['files'].get('aln_fasta')
    aln_sto_intgz = contents['files'].get('aln_sto')

    if aln_fasta_intgz and aln_sto_intgz:
        # Both version of the alignment are in the refpkg
        with open('refpkg.aln.fasta','w') as out_aln_fasta_h:
            out_aln_fasta_h.write(
                tar_h.extractfile(
                    tar_contents_dict[aln_fasta_intgz]
                ).read().decode('utf-8')
            )
        with open('refpkg.aln.sto','w') as out_aln_sto_h:
            out_aln_sto_h.write(
                tar_h.extractfile(
                    tar_contents_dict[aln_sto_intgz]
                ).read().decode('utf-8')
            )
    elif aln_fasta_intgz:
        # Only fasta exists
        with open('refpkg.aln.fasta','w') as out_aln_fasta_h:
            out_aln_fasta_h.write(
                tar_h.extractfile(
                    tar_contents_dict[aln_fasta_intgz]
                ).read().decode('utf-8')
            )
        # And convert to sto format
        with open('refpkg.aln.sto','w') as out_aln_sto_h:
            AlignIO.write(
                AlignIO.read(
                    tar_h.extractfile(tar_contents_dict[aln_fasta_intgz]),
                    'fasta'),
                out_aln_sto_h,
                'stockholm'
            )
    elif aln_sto_intgz:
        # Only STO exists
        with open('refpkg.aln.sto','w') as out_aln_sto_h:
            out_aln_sto_h.write(
                tar_h.extractfile(
                    tar_contents_dict[aln_sto_intgz]
                ).read().decode('utf-8')
            )
        with sopen('refpkg.aln.fasta','w') as out_aln_fasta_h:
            AlignIO.write(
                AlignIO.read(
                                tar_h.extractfile(tar_contents_dict[aln_sto_intgz]),
                                'stockholm'),
                out_aln_fasta_h,
                'fasta'
            )
    else:
        # NO alignment present
        raise Exception("Refset does not contain an alignment")
    """
}

process combineAln_SV_refpkg {
    container = 'golob/infernal:1.1.2_bcw_0.2.0'
    label = 'mem_veryhigh'

    input:
        file sv_aln_sto_f 
        file refpkg_aln_sto_f
        
    
    output:
        file "sv_refpkg.aln.sto" into sv_refpkg_aln_sto_f
    
    """
    esl-alimerge --dna \
     -o sv_refpkg.aln.sto \
     ${sv_aln_sto_f} ${refpkg_aln_sto_f}
    """
}


//  Step 3.c. Place SV via pplacer
process pplacerPlacement {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'mem_veryhigh'

    publishDir "${params.output}/placement", mode: 'copy'

    input:
        file sv_refpkg_aln_sto_f
        file refpkg_tgz_f
    output:
        file 'dedup.jplace' into dedup_jplace_f
    
    afterScript "rm -rf refpkg/"
    """
    mkdir -p refpkg/ &&
    tar xzvf ${refpkg_tgz_f} -C ./refpkg &&
    pplacer -p -j ${task.cpus} \
    --inform-prior --prior-lower ${params.pplacer_prior_lower} --map-identity \
    -c refpkg/ ${sv_refpkg_aln_sto_f} \
    -o dedup.jplace
    """
}

//  Step 3.d. Reduplicate placements
process pplacerReduplicate {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'

    publishDir "${params.output}/placement", mode: 'copy'

    input:
        file dedup_jplace_f
        file sv_weights_f
    output:
        file 'redup.jplace.gz'
    
    """
    guppy redup -m \
    -o /dev/stdout \
    -d ${sv_weights_f} \
    ${dedup_jplace_f} \
    | gzip > redup.jplace.gz
    """
}


//  Step 3.e. ADCL metric
process pplacerADCL {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'

    publishDir "${params.output}/placement", mode: 'copy'

    input:
        file dedup_jplace_f
    output:
        file 'adcl.csv.gz'
    
    """
    (echo name,adcl,weight && 
    guppy adcl --no-collapse ${dedup_jplace_f} -o /dev/stdout) | 
    gzip > adcl.csv.gz
    """
}

//  Step 3.f. EDPL metric
process pplacerEDPL {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'

    publishDir "${params.output}/placement", mode: 'copy'

    input:
        file dedup_jplace_f
    output:
        file 'edpl.csv.gz'
    
    """
    (echo name,edpl && guppy edpl --csv ${dedup_jplace_f} -o /dev/stdout) | 
    gzip > edpl.csv.gz
    """
}

//  Step 3.g. (e/l)PCA
process pplacerPCA {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'
    afterScript "rm -r refpkg/"
    publishDir "${params.output}/placement", mode: 'copy'
    errorStrategy = 'ignore'

    input:
        file refpkg_tgz_f
        file dedup_jplace_f
        file sv_map_f
    output:
        file 'pca/epca.proj'
        file 'pca/epca.xml'
        file 'pca/epca.trans'
        file 'pca/lpca.proj'
        file 'pca/lpca.xml'
        file 'pca/lpca.trans'
    
    """
    mkdir -p refpkg/ && mkdir -p pca/
    tar xzvf ${refpkg_tgz_f} -C refpkg/ &&
    guppy epca ${dedup_jplace_f}:${sv_map_f} -c refpkg/ --out-dir pca/ --prefix epca &&
    guppy lpca ${dedup_jplace_f}:${sv_map_f} -c refpkg/ --out-dir pca/ --prefix lpca
    """
}

//  Step 3.h. Alpha diversity
process pplacerAlphaDiversity {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'

    publishDir "${params.output}/placement", mode: 'copy'

    input:
        file dedup_jplace_f
        file sv_map_f
    output:
        file 'alpha_diversity.csv.gz'

    
    """
    guppy fpd --csv --include-pendant --chao-d 0,1,1.00001,2,3,4,5 \
    ${dedup_jplace_f}:${sv_map_f} |
    gzip > alpha_diversity.csv.gz
    """
}

//  Step 3.i. KR (phylogenetic) distance 
process pplacerKR {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'
    afterScript "rm -r refpkg/"
    publishDir "${params.output}/placement", mode: 'copy'

    input:
        file refpkg_tgz_f
        file dedup_jplace_f
        file sv_map_f
    output:
        file 'kr_distance.csv.gz'

    
    """
    mkdir -p refpkg/
    tar xzvf ${refpkg_tgz_f} -C refpkg/
    guppy kr --list-out -c refpkg/ ${dedup_jplace_f}:${sv_map_f} |
    gzip > kr_distance.csv.gz
    """
}
//
//  END STEP 3: Placement
//

//
//  START STEP 4: Classification
//

params.pp_classifer = 'hybrid2'

params.pp_likelihood_cutoff = 0.9
params.pp_bayes_cutoff = 1.0
params.pp_multiclass_min = 0.2
params.pp_bootstrap_cutoff = 0.8
params.pp_bootstrap_extension_cutoff = 0.4
params.pp_nbc_boot = 100
params.pp_nbc_target_rank = 'genus'
params.pp_nbc_word_length = 8
params.pp_seed = 1

// Step 4.a. Prep the placement DB
process classifyDB_Prep {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'
    afterScript "rm -r refpkg/"
    cache = false

    input:
        file refpkg_tgz_f
        file sv_map_f
    
    output:
        file 'classify.prep.db' into classify_db_prepped
    

    """
    mkdir -p refpkg/
    tar xzvf ${refpkg_tgz_f} -C refpkg/
    rppr prep_db -c refpkg/ --sqlite classify.prep.db
    (echo "name,specimen"; cat ${sv_map_f}) |
    csvsql --table seq_info --insert --snifflimit 1000 --db sqlite:///classify.prep.db
    """
}

// Step 4.b. Classify SV
process classifySV {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'mem_veryhigh'
    afterScript "rm -r refpkg/"
    cache = false

    input:
        file refpkg_tgz_f
        file classify_db_prepped
        file dedup_jplace_f
        file sv_refpkg_aln_sto_f
    
    output:
        file 'classify.classified.db' into classifyDB_classified

    """
    mkdir -p refpkg/
    tar xzvf ${refpkg_tgz_f} -C refpkg/
    guppy classify --pp \
    --classifier ${params.pp_classifer} \
    -j ${task.cpus} \
    -c refpkg/ \
    --nbc-sequences ${sv_refpkg_aln_sto_f} \
    --sqlite ${classify_db_prepped} \
    --seed ${params.pp_seed} \
    --cutoff ${params.pp_likelihood_cutoff} \
    --bayes-cutoff ${params.pp_bayes_cutoff} \
    --multiclass-min ${params.pp_multiclass_min} \
    --bootstrap-cutoff ${params.pp_bootstrap_cutoff} \
    --bootstrap-extension-cutoff ${params.pp_bootstrap_extension_cutoff} \
    --word-length ${params.pp_nbc_word_length} \
    --nbc-rank ${params.pp_nbc_target_rank} \
    --n-boot ${params.pp_nbc_boot} \
    ${dedup_jplace_f}
    cp ${classify_db_prepped} classify.classified.db
    """
}

// Step 4.c. Concatenate placements
process classifyMCC {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'
    cache = false
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        file classifyDB_classified
        file sv_weights_f

    output:
        file 'classify.mcc.db' into classifyDB_mcc

    """
    multiclass_concat.py -k \
    --dedup-info ${sv_weights_f} ${classifyDB_classified}
    cp ${classifyDB_classified} classify.mcc.db
    """
}

// Step 4.d. Tabular outputs

Channel.from(
    'phylum', 'class', 'order', 'family', 'genus', 'species'
).set { classify_ranks }

classify_ranks.combine(
    classifyDB_mcc
).combine(
    sv_map_f
).set { classify_rank_ch }

process classifyTables {
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'
    label = 'io_limited'
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        set val(rank), file(classifyDB_mcc), file(sv_map_for_tables_f) from classify_rank_ch

    output:
        set val(rank), file("tables/by_specimen.${rank}.csv"), file("tables/by_taxon.${rank}.csv"), file("tables/tallies_wide.${rank}.csv") into classify_rank_tables_ch

    """
    mkdir -p tables/
    python2 /usr/bin/classif_table.py ${classifyDB_mcc} \
    tables/by_taxon.${rank}.csv \
    --rank ${rank} \
    --specimen-map ${sv_map_for_tables_f} \
    --by-specimen tables/by_specimen.${rank}.csv \
    --tallies-wide tables/tallies_wide.${rank}.csv
    """
}


//
//  END STEP 4: Classification
//


// */

