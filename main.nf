#!/usr/bin/env nextflow

/*
    Maliampi via nextflow!
    Steps in a typical 16S rRNA run:
    1) Make sequence variants (with dada2)
    2) create (vs load) a reference package
    3) Place SV on the reference package
    4) Classify the SV using the placements + reference package
*/
nextflow.preview.dsl=2

// Defaults for parameters
params.help = false
params.sv_only = false
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
params.errM_randomize = 'TRUE'
params.errM_nbases = '1e8'
params.chimera_method = 'consensus'

// Good's filtering
params.goods_convergence = 0.0001
params.min_sv_prev = 2
params.goods_min_reads = 30

// Refpkg
params.repo_min_id = 0.8
params.repo_max_accepts = 10
params.cmalign_mxsize = 8196
params.raxml_model = 'GTRGAMMA'
params.raxml_parsiomony_seed = 12345
params.taxdmp = false





// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/maliampi <ARGUMENTS>
    
    Required Arguments:
        --manifest            CSV file listing samples
                                At a minimum must have columns:
                                    specimen: A unique identifier 
                                    R1: forward read
                                    R2: reverse read fq

                                optional columns:
                                    batch: sequencing / library batch. Should be filename safe
                                    I1: forward index file (for checking demultiplexing)
                                    I2: reverse index file
        --repo_fasta          Repository of 16S rRNA genes.
        --repo_si             Information about the 16S rRNA genes.
        --email               Email (for NCBI)
    Options:
      Common to all:
        --output              Directory to place outputs (default invocation dir)
                                Maliampi will create a directory structure under this directory
        -w                    Working directory. Defaults to `./work`
        -resume                 Attempt to restart from a prior run, only completely changed steps

    Flow options:
        --sv_only             Stop after making sequence variants   

    SV-DADA2 options:
        --trimLeft              How far to trim on the left (default = 0)
        --maxN                  (default = 0)
        --maxEE                 (default = Inf)
        --truncLenF             (default = 0)
        --truncLenR             (default = 0)
        --truncQ                (default = 2)
    
    Good's Coverage Filtering
        --goods_convergence     Change in Good's converage to consider converged (default = 0.0001)
        --min_sv_prev           Minimum prevalance for an SV to be kept (default = 2)
        --goods_min_reads       Minimum reads to conisder before Good's is considered converged

    Ref Package required:
        --repo_fasta            FASTA file containing reference sequences (required)
        --repo_si               CSV file with information about the repo reads (required)
    Ref Package options (defaults generally fine):
        --repo_min_id           Minimum percent ID to a SV to be recruited (default = 0.8)
        --repo_max_accepts      Maximum number of recruits per SV (default = 10)
        --cmalign_mxsize        Infernal cmalign mxsize (default = 8196)
        --raxml_model           RAxML model for tree formation (default = 'GTRGAMMA')
        --raxml_parsiomony_seed (default = 12345)
        --taxdmp                Path to taxdmp.zip. If not provided, it will be downloaded

    Placement / Classification Options (defaults generally fine):
        --pp_classifer                  pplacer classifer (default = 'hybrid2')
        --pp_likelihood_cutoff          (default = 0.9)
        --pp_bayes_cutoff               (default = 1.0)
        --pp_multiclass_min             (default = 0.2)
        --pp_bootstrap_cutoff           (default = 0.8)
        --pp_bootstrap_extension_cutoff (default = 0.4)
        --pp_nbc_boot                   (default = 100)
        --pp_nbc_target_rank            (default = 'genus')
        --pp_nbc_word_length            (default = 8)
        --pp_seed                       (default = 1)
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.manifest == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Modules
include read_manifest from './modules/manifest'
include output_failed from './modules/preprocess' params (
    output: params.output
)
include preprocess_wf from './modules/preprocess'
include dada2_wf from './modules/dada2' params (
    output: params.output,
    trimLeft: params.trimLeft,
    maxN: params.maxN,
    maxEE: params.maxEE,
    truncLenF: params.truncLenF,
    truncLenR: params.truncLenR,
    truncQ: params.truncQ,
    errM_maxConsist: params.errM_maxConsist,
    errM_randomize: params.errM_randomize,
    errM_nbases: params.errM_nbases,
    chimera_method: params.chimera_method,
    goods_convergence: params.goods_convergence,
    min_sv_prev: params.min_sv_prev,
    goods_min_reads :params.goods_min_reads
)

include make_refpkg_wf from './modules/refpackage' params (
    output: params.output,
    repo_fasta: params.repo_fasta,
    repo_si: params.repo_si,
    repo_min_id: params.repo_min_id,
    repo_max_accepts: params.repo_max_accepts,
    cmalign_mxsize: params.cmalign_mxsize,
    raxml_model: params.raxml_model,
    raxml_parsiomony_seed: params.raxml_parsiomony_seed,
    taxdmp: params.taxdmp

)

// STEP 0: Read manifest and verify files.

workflow {
    main:

    //
    //  Step 0: Load manifest and preprocess
    //

    // Load manifest!
    manifest = read_manifest(
        Channel.from(
            file(params.manifest)
        )
    )
    // manifest.valid_paired_indexed contains indexed paired reads
    // manifest.valid_paired contains pairs verified to exist but without index.

    // Preprocess
    preprocess_wf(
        manifest.valid_paired_indexed,
        manifest.valid_paired
    )        
    // preprocess_wf.out.valid is the reads that survived the preprocessing steps.
    // preprocess_wf.out.empty are the reads that ended up empty with preprocessing

    //
    // Step 1: DADA2 to make sequence variants.
    //

    dada2_wf(preprocess_wf.out.valid)

    //
    // Report specimens that failed at any step of making SVs
    //

    output_failed(            
        manifest.other.map { [it.specimen, 'failed at manifest'] }.mix(
        preprocess_wf.out.empty.map{ [it[0], 'preprocessing'] }).mix(
        dada2_wf.out.failures)
        .toList()
        .transpose()
        .toList()
    )
   
    if (params.sv_only) {
        return;
    }

    // implicit else

    //
    //  STEP 2: Reference package
    //

    make_refpkg_wf(
        dada2_wf.out.sv_fasta
    )

    
}
return


if (params.sv_only == false) {



//
//  START STEP 3: Placement
//
params.pplacer_prior_lower = 0.01

//  Step 3.a. Align SV

process alignSV {
    container = "${container__infernal}"
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
    container = "${container__fastatools}"
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
    container = "${container__infernal}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
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
    container = "${container__pplacer}"
    label = 'io_limited'
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        set val(rank), file(classifyDB_mcc), file(sv_map_for_tables_f) from classify_rank_ch

    output:
        set val(rank), file("tables/by_specimen.${rank}.csv"), file("tables/by_taxon.${rank}.csv"), file("tables/tallies_wide.${rank}.csv") into classify_rank_tables_ch

    """
    mkdir -p tables/
    classif_table.py ${classifyDB_mcc} \
    tables/by_taxon.${rank}.csv \
    --rank ${rank} \
    --specimen-map ${sv_map_for_tables_f} \
    --by-specimen tables/by_specimen.${rank}.csv \
    --tallies-wide tables/tallies_wide.${rank}.csv
    """
}

} // END NOT SV ONLY

//
//  END STEP 4: Classification
//


// */

