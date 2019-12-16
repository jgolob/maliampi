#!/usr/bin/env nextflow

/*
    Maliampi via nextflow!
    Steps in a typical 16S rRNA run:
    1) Make sequence variants (with dada2)
    2) create (vs load) a reference package
    3) Place SV on the reference package
    4) Classify the SV using the placements + reference package

    This subpackage only does 2, 3 and 4.
*/


// container versions!
container__barcodecop = "golob/barcodecop:0.4.1__bcw_0.3.0"
container__dada2 = "golob/dada2:1.12.0.ub.1804__bcw.0.3.1"
container__fastcombineseqtab = "golob/dada2-fast-combineseqtab:0.5.0__1.12.0__BCW_0.3.1"
container__dada2pplacer = "golob/dada2-pplacer:0.8.0__bcw_0.3.1A"
container__vsearch = "golob/vsearch:2.7.1_bcw_0.2.0"
container__fastatools = "golob/fastatools:0.7.1__bcw.0.3.1"
container__pplacer = "golob/pplacer:1.1alpha19rc_BCW_0.3.1A"
container__seqinfosync = "golob/seqinfo_taxonomy_sync:0.2.1__bcw.0.3.0"
container__infernal = "golob/infernal:1.1.2_bcw_0.3.1"
container__raxml = "golob/raxml:8.2.11_bcw_0.3.0"

// Defaults for parameters
params.help = false


// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/maliampi <ARGUMENTS>
    
    Required Arguments:
        --fasta               Sequence Variant FASTA
        --weights             CSV file with weights by SV
        --map                 CSV file with map of SV <-> community
        --repo_fasta          Repository FASTA
        --repo_si             Repository Seq Info (CSV)
    Options:
      Common to all:
        --output              Directory to place outputs (default invocation dir)
                                Maliampi will create a directory structure under this directory
        -w                    Working directory. Defaults to `./work`
        -resume                 Attempt to restart from a prior run, only completely changed steps


    Placement / Classification Options (defaults generally fine):
        --cmalign_mxsize                Infernal cmalign mxsize (default = 8196)
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
if (params.fasta == null || params.weights == null || params.map == null || params.repo_fasta == null || params.repo_si == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}


// Load inputs

sv_fasta_f = Channel.value(file(params.fasta))
sv_weights_f = Channel.value(file(params.weights))
sv_map_f = Channel.value(file(params.map))

//
//  START STEP 2: Reference package
//
params.repo_min_id = 0.8
params.repo_max_accepts = 10
params.cmalign_mxsize = 8196
params.raxml_model = 'GTRGAMMA'
params.raxml_parsiomony_seed = 12345
params.taxdmp = false


// Step 2.a. Search the repo for matches
// load the repo
Channel.value(file(params.repo_fasta))
    .set{ repo_fasta}
Channel.value(file(params.repo_si))
    .set{ repo_si }

process refpkgSearchRepo {
    container "${container__vsearch}"
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
    container = "${container__fastatools}"
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

if ( (params.taxdmp == false) || file(params.taxdmp).isEmpty() ) {
    process DlBuildTaxtasticDB {
        container = "${container__pplacer}"
        label = 'io_limited'
        // errorStrategy = 'retry'

        output:
            file "taxonomy.db" into taxonomy_db_f

        afterScript "rm -rf dl/"


        """
        mkdir -p dl/ && \
        taxit new_database taxonomy.db -u ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip -p dl/
        """

    }

} else {
    taxdump_zip_f = file(params.taxdmp)
    process buildTaxtasticDB {
        container = "${container__pplacer}"
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
    container = "${container__seqinfosync}"
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
    container = "${container__infernal}"
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
    container = "${container__fastatools}"
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
    container = "${container__raxml}"
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

// Step 2.xx Remove cruft from tree stats

process raxmlTree_cleanupInfo {
    container = "${container__raxml}"
    label = 'io_limited'
    errorStrategy = 'retry'

    input:
        file "RAxML_info.unclean.refpkg" from refpkg_tree_stats_f
    
    output:
        file "RAxML_info.refpkg" into refpkg_tree_stats_clean_f


"""
#!/usr/bin/env python
with open("RAxML_info.refpkg",'wt') as out_h:
    with open("RAxML_info.unclean.refpkg", 'rt') as in_h:
        past_cruft = False
        for l in in_h:
            if "This is RAxML version" == l[0:21]:
                past_cruft = True
            if past_cruft:
                out_h.write(l)
"""
}

// Step 2.xx Make a tax table for the refpkg sequences
process taxtableForSI {
    container = "${container__pplacer}"
    label = 'io_limited'
    // errorStrategy = 'retry'

    input:
        file taxonomy_db_f 
        file refpkg_si_corr_f
    output:
        file "refpkg.taxtable.csv" into refpkg_tt_f

    """
    taxit taxtable ${taxonomy_db_f} \
    --seq-info ${refpkg_si_corr_f} \
    --outfile refpkg.taxtable.csv
    """
}

// Step 2.xx Obtain the CM used for the alignment
process obtainCM {
    container = "${container__infernal}"
    label = 'io_limited'

    output:
        file "refpkg.cm" into refpkg_cm
    
    """
    cp /cmalign/data/SSU_rRNA_bacteria.cm refpkg.cm
    """
}

// Step 2.xx Combine into a refpkg
process combineRefpkg {
    container = "${container__pplacer}"
    label = 'io_limited'

    afterScript("rm -rf refpkg/*")
    publishDir "${params.output}/refpkg/", mode: 'copy'

    input:
        file recruits_aln_fasta_f
        file recruits_aln_sto_f
        file refpkg_tree_f 
        file refpkg_tree_stats_clean_f 
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
    --tree-stats ${refpkg_tree_stats_clean_f} \
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

//
//  END STEP 4: Classification
//


// */

