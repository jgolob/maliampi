#!/usr/bin/env nextflow

/*
    Maliampi via nextflow!
    Steps in a typical 16S rRNA run:
    1) Make sequence variants (with dada2)
    2) create (vs load) a reference package
    3) Place SV on the reference package
    4) Derive taxonomy, phylotypes, and placement statistics
*/
nextflow.enable.dsl=2

// Defaults for parameters
params.help = false
params.sv_only = false
params.skip_taxonomy = false
params.skip_stats = false
params.skip_phylotypes = false
params.placer = 'epang'
params.legacy_redup = false
params.refpkg = null
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
        --project_id          Stable project identity for every observation
        --dataset_id          Stable dataset identity for this SV artifact
        One reference source:
          --refpkg            Existing reference package, or
          --repo_fasta        Repository of 16S rRNA genes
          --repo_si           Information about the reference sequences
          --email             Contact address for refpkg construction
    Options:
      Common to all:
        --output              Directory to place outputs (default invocation dir)
                                Maliampi will create a directory structure under this directory
        -w                    Working directory. Defaults to `./work`
        -resume                 Attempt to restart from a prior run, only completely changed steps
        --sv_only              Stop after canonical SV H5AD generation
        --skip_taxonomy        Do not generate taxonomy artifacts
        --skip_phylotypes      Do not generate phylotype artifacts
        --skip_stats           Do not generate placement statistics

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

    Placement Options:
        --placer                        Placement engine: epang (default) or pplacer
        --legacy_redup                  Emit legacy pplacer-reduplicated output
    """.stripIndent()
}

// Modules
include { make_refpkg_wf } from './modules/refpackage'
include { sv } from './subworkflows/local/sv'
include { place } from './subworkflows/local/place'
include { taxonomy } from './subworkflows/local/taxonomy'
include { stats } from './subworkflows/local/stats'
include { phylotypes } from './subworkflows/local/phylotypes'
include { effective_phylotype_thresholds } from './subworkflows/local/phylotypes'
include { WriteSvRegistry } from './modules/sv_h5ad'

// STEP 0: Read manifest and verify files.

workflow {
    main:
    // Show help message if the user specifies the --help flag at runtime
    if (
        params.help || 
        (params.manifest == null) ||
        (params.project_id == null) ||
        (params.dataset_id == null) ||
        (!params.sv_only && params.refpkg == null && (params.repo_fasta == null || params.repo_si == null || params.email == null))
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

    sv(channel.fromPath(params.manifest, checkIfExists: true))

    if (!params.sv_only) {
        //
        //  STEP 2: Reference package
        //

        if (params.refpkg != null) {
            refpkg = channel.fromPath(params.refpkg, checkIfExists: true)
        } else {
            WriteSvRegistry(sv.out.sv_h5ad)
            make_refpkg_wf(sv.out.sv_fasta, WriteSvRegistry.out.registry)
            refpkg = make_refpkg_wf.out.refpkg_tgz
        }

    //
    // STEP 3. Place and Classify
    //
        place(
            sv.out.sv_h5ad,
            refpkg,
            params.placer,
            params.legacy_redup,
        )

        if (!params.skip_taxonomy) {
        // Taxonomy validates the same immutable source package independently.
        // Do not feed it placement's multi-output validation process channel.
            taxonomy(place.out.dedup_jplace, refpkg, sv.out.sv_h5ad, place.out.placement_identity)
        }
        if (!params.skip_stats) {
            stats(place.out.dedup_jplace, place.out.specimen_jplaces)
        }
        if (!params.skip_phylotypes) {
        def thresholds = effective_phylotype_thresholds(
            params.phylotype_thresholds,
            params.phylotype_add_thresholds,
            params.phylotype_distance,
            params.phylotype_kr_thresholds,
        )
            phylotypes(
                place.out.dedup_jplace,
                sv.out.sv_h5ad,
                place.out.placement_identity,
                params.dataset_id,
                thresholds,
                params.phylotype_distance,
                params.phylotype_lwr_overlap,
            )
        }
    }
    
}
