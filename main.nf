#!/usr/bin/env nextflow

/*
    Maliampi via nextflow!
    Steps in a typical 16S rRNA run:
    1) Make sequence variants (with dada2)
    2) create (vs load) a reference package
    3) Place SV on the reference package
    4) Classify the SV using the placements + reference package
*/
nextflow.enable.dsl=2

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
params.raxmlng_model = 'GTR+G'
params.raxmlng_parsimony_trees = 1
params.raxmlng_random_trees = 1
params.raxmlng_bootstrap_cutoff = 0.3
params.raxmlng_seed = 12345
params.raxml_model = 'GTRGAMMA'
params.raxml_parsiomony_seed = 12345
params.raxml = 'og'
params.taxdmp = false

// pplacer place
params.pplacer_prior_lower = 0.01

// pplacer classify
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

    SV-DADA2 options:
        --trimLeft              How far to trim on the left (default = 0)
        --maxN                  (default = 0)
        --maxEE                 (default = Inf)
        --truncLenF             (default = 0)
        --truncLenR             (default = 0)
        --truncQ                (default = 2)
        --minOverlap            (default = 12)
        --maxMismatch           (default = 0)


    Ref Package options (defaults generally fine):
        --raxml                     Which raxml to use: og (original) or ng (new). Default: og
        --repo_min_id               Minimum percent ID to a SV to be recruited (default = 0.8)
        --repo_max_accepts          Maximum number of recruits per SV (default = 10)
        --cmalign_mxsize            Infernal cmalign mxsize (default = 8196)
        --raxml_model               RAxML model for tree formation (default = 'GTRGAMMA')
        --raxml_parsiomony_seed     (default = 12345)        
        --raxmlng_model             Subsitution model (default 'GTR+G')
        --raxmlng_parsimony_trees   How many seed parsimony trees (default 1)
        --raxmlng_random_trees      How many seed random trees (default 1)
        --raxmlng_bootstrap_cutoff  When to stop boostraps (default = 0.3)
        --raxmlng_seed              Random seed for RAxML-ng (default = 12345)
        --taxdmp                    (Optional) taxdmp.zip from the repository

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

// Modules
include { read_manifest } from './modules/manifest'
include { output_failed } from './modules/preprocess'
include { preprocess_wf } from './modules/preprocess'
include { dada2_wf } from './modules/dada2'
include { make_refpkg_wf } from './modules/refpackage'
include { epang_place_classify_wf } from './modules/epang_place_classify'

// STEP 0: Read manifest and verify files.

workflow {
    main:
    // Show help message if the user specifies the --help flag at runtime
    if (
        params.help || 
        (params.manifest == null) ||
        (params.repo_fasta == null) ||
        (params.repo_si == null) ||
        (params.email == null)
    ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    } 
    // Implicit else

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
        manifest.valid_paired,
        manifest.valid_unpaired
    )        
    // preprocess_wf.out.valid is the reads that survived the preprocessing steps.
    // preprocess_wf.out.empty are the reads that ended up empty with preprocessing

    //
    // Step 1: DADA2 to make sequence variants.
    //

    dada2_wf(
        preprocess_wf.out.miseq_pe,
        preprocess_wf.out.miseq_se,
        preprocess_wf.out.pyro
    )

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

    //
    //  STEP 2: Reference package
    //

    make_refpkg_wf(
        dada2_wf.out.sv_fasta
    )

    //
    // STEP 3. Place and Classify
    //
    epang_place_classify_wf(
        dada2_wf.out.sv_fasta,
        make_refpkg_wf.out.refpkg_tgz,
        dada2_wf.out.sv_long
    )
    
}

