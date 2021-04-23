//
//  PPlacer Place and Classify
//
nextflow.preview.dsl=2

// Paramteters
// common
params.output = '.'
params.help = false

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
params.cmalign_mxsize = 8196

// Containers!
container__infernal = "golob/infernal:1.1.2_bcw_0.3.1"
container__fastatools = "golob/fastatools:0.7.1__bcw.0.3.1"
container__pplacer = "golob/pplacer:1.1alpha19rc_BCW_0.3.1A"
container__dada2pplacer = "golob/dada2-pplacer:0.8.0__bcw_0.3.1A"
container__epang = "quay.io/biocontainers/epa-ng:0.3.8--h9a82719_1"
container__gappa = 'quay.io/biocontainers/gappa:0.7.1--h9a82719_1'


workflow place_classify_wf {
    take:
        sv_fasta_f
        refpkg_tgz_f
        sv_weights_f
        sv_map_f
        sv_long_f

    main:

    //
    // Step 1. Align the SV
    // 
    AlignSV(
        sv_fasta_f
    )

    //
    //  PLACEMENT
    //

    //
    //  Step 2. Extract reference alignment, and combine
    //
    
    // Step 2a. Extract bits from the reference
    ExtractRefpkgForEPAng(
        refpkg_tgz_f
    )

    // Step 2b. Combine SV and refpkg alignments
    CombineAln_SV_refpkg(
        AlignSV.out[0],
        ExtractRefpkgForEPAng.out.ref_aln_sto
    )

    // Step 2c. Convert combined alignment to FASTA format
    ConvertAlnToFasta(
        CombineAln_SV_refpkg.out
    )

    //
    //  Step 3. Place SV via epa-ng
    //

    EPAngPlacement(
        ExtractRefpkgForEPAng.out.ref_aln_fasta,
        ConvertAlnToFasta.out,
        ExtractRefpkgForEPAng.out.model,
        ExtractRefpkgForEPAng.out.tree
    )

    //
    //  Step 4. Reduplicate placements
    //

    MakeSplit(
        sv_long_f
    )

    GappaSplit(
        EPAngPlacement.out,
        MakeSplit.out
    )

    //
    //  Step 5. ADCL metric
    //
    PplacerADCL(
        EPAngPlacement.out
    )

    //
    //  Step 6. EDPL metric
    //
    EDPL(
        EPAngPlacement.out
    )

    //
    //  Step 7. xPCA
    //

    Gappa_ePCA(
         GappaSplit.out
    )

    //
    //  Step 8. Alpha diversity
    //
    PplacerAlphaDiversity(
        EPAngPlacement.out,
        sv_map_f 
    )

    //
    //  Step 9. KR (phylogenetic) distance 
    //
    Gappa_KRD(
        GappaSplit.out
    )
    //
    //  END Placement
    //

    //
    //  CLASSIFY
    //

    MakeEPAngTaxonomy(
        ExtractRefpkgForEPAng.out.leaf_info,
        ExtractRefpkgForEPAng.out.taxonomy,
    )

    Gappa_Classify(
        EPAngPlacement.out,
        MakeEPAngTaxonomy.out
    )

    emit:
        jplace_dedup = EPAngPlacement.out
        //jplace_redup = PplacerReduplicate.out

}

process AlignSV {
    container = "${container__infernal}"
    label = 'mem_veryhigh'

    input:
        file sv_fasta_f
    
    output:
        file "sv.aln.sto"
        file "sv.aln.scores"
        
    
    """
    cmalign \
    --cpu ${task.cpus} --noprob --dnaout --mxsize ${params.cmalign_mxsize} \
    --sfile sv.aln.scores -o sv.aln.sto \
    /cmalign/data/SSU_rRNA_bacteria.cm ${sv_fasta_f}
    """
}

process CombineAln_SV_refpkg {
    container = "${container__infernal}"
    label = 'mem_veryhigh'

    input:
        file sv_aln_sto_f 
        file refpkg_aln_sto_f
        
    
    output:
        file "sv_refpkg.aln.sto"
    
    """
    esl-alimerge --dna \
     -o sv_refpkg.aln.sto \
     ${sv_aln_sto_f} ${refpkg_aln_sto_f}
    """
}

process ConvertAlnToFasta {
    container = "${container__fastatools}"
    label = 'io_limited'
    errorStrategy "retry"

    input: 
        file combined_aln_sto_f
    
    output:
        file "combined.aln.fasta"
    
    """
    #!/usr/bin/env python
    from Bio import AlignIO

    with open('combined.aln.fasta', 'wt') as out_h:
        AlignIO.write(
            AlignIO.read(
                open('${combined_aln_sto_f}', 'rt'),
                'stockholm'
            ),
            out_h,
            'fasta'
        )
    """
}

process ExtractRefpkgForEPAng {
    container = "${container__fastatools}"
    label = 'io_limited'
    
    input:
        file refpkg_tgz_f

    output:
        path 'refpkg_tree.nwk', emit: tree
        path 'refpkg.aln.fasta', emit: ref_aln_fasta
        path 'refpkg.aln.sto', emit: ref_aln_sto
        path 'model.txt', emit: model
        path 'leaf_info.csv', emit: leaf_info
        path 'taxonomy.csv', emit: taxonomy

"""
#!/usr/bin/env python
import tarfile
import json
import os
import re

tar_h = tarfile.open('${refpkg_tgz_f}')
tar_contents_dict = {os.path.basename(f.name): f for f in tar_h.getmembers()}
contents = json.loads(
    tar_h.extractfile(
        tar_contents_dict['CONTENTS.json']
    ).read().decode('utf-8')
)

with open('refpkg_tree.nwk', 'wt') as tree_h:
    tree_h.writelines(
        tar_h.extractfile(
                tar_contents_dict[contents['files'].get('tree')]
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
# Model
phylo_model = json.loads(tar_h.extractfile(
        tar_contents_dict[contents['files'].get('phylo_model')]
    ).read().decode('utf-8')
).get('subs_rates')

re_basefreq = re.compile(r'Base frequencies: (?P<A>0\\.\\d+) (?P<C>0\\.\\d+) (?P<G>0\\.\\d+) (?P<T>0\\.\\d+)')
bf_m = re_basefreq.search(tar_h.extractfile(
        tar_contents_dict[contents['files'].get('tree_stats')]
    ).read().decode('utf-8'))
with open('model.txt', 'wt') as model_h:
    model_h.writelines( 
        "GTR{"+
        "{}/{}/{}/{}/{}/{}".format(
            phylo_model['ac'],
            phylo_model['ag'],
            phylo_model['at'],
            phylo_model['cg'],
            phylo_model['ct'],
            phylo_model['gt'],
        )
        +"}"+"+FU{"+
        "{}/{}/{}/{}".format(
            bf_m['A'],
            bf_m['C'],
            bf_m['G'],
            bf_m['T'],
        )
        +"}"
    )

with open('leaf_info.csv', 'wt') as leaf_h:
    leaf_h.write(tar_h.extractfile(
        tar_contents_dict[contents['files'].get('seq_info')]
    ).read().decode('utf-8'))

with open('taxonomy.csv', 'wt') as leaf_h:
    leaf_h.write(tar_h.extractfile(
        tar_contents_dict[contents['files'].get('taxonomy')]
    ).read().decode('utf-8'))


"""

}

process EPAngPlacement {
    container = "${container__epang}"
    label = 'mem_veryhigh'
    publishDir "${params.output}/placement", mode: 'copy'
    input:
        file refpkg_aln_fasta
        file combined_aln_fasta
        file model
        file ref_tree

    output:
        file 'dedup.jplace'
    """
    set -e

    epa-ng --split ${refpkg_aln_fasta} ${combined_aln_fasta}
    model=`cat ${model}`
    
    epa-ng -t ${ref_tree} \
    -s reference.fasta -q query.fasta \
    -m \$model -T ${task.cpus} \
    --baseball-heur

    mv epa_result.jplace dedup.jplace
    """
}


process MakeEPAngTaxonomy {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/refpkg", mode: 'copy'

    input:
        path leaf_info_f
        path taxonomy_f

    output:
        path 'epang_taxon_file.tsv'

"""
#!/usr/bin/env python
import csv

tax_dict = {
    r['tax_id']: r for r in
    csv.DictReader(open('${taxonomy_f}', 'rt'))
}
tax_names = {
    tax_id: r['tax_name']
    for tax_id, r in tax_dict.items()
}
RANKS = [
    'superkingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species',
]
with open('epang_taxon_file.tsv', 'wt') as tf_h:
    tf_w = csv.writer(tf_h, delimiter='\\t')
    for row in csv.DictReader(open('${leaf_info_f}', 'rt')):
        tax_id = row.get('tax_id', None)
        if tax_id is None:
            continue
        # Implicit else
        tax_lineage = tax_dict.get(tax_id, None)
        if tax_lineage is None:
            continue
        lineage_str = ";".join([
            tax_names.get(tax_lineage.get(rank, ""), "")
            for rank in RANKS
        ])
        tf_w.writerow([row['seqname'], lineage_str])

"""

}

process Gappa_Classify {
    container = "${container__gappa}"
    label = 'multithreaded'
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        path dedup_jplace
        path taxon_file
    
    output:
        path 'per_query.tsv'

    """
    set -e


    gappa examine assign \
    --per-query-results \
    --verbose \
    --threads ${task.cpus} \
    --jplace-path ${dedup_jplace} \
    --taxon-file ${taxon_file} \
    
    """

}

process EDPL {
    container = "${container__gappa}"
    label = 'multithreaded'
    publishDir "${params.output}/placement", mode: 'copy'

    input:
        path dedup_jplace
    
    output:
        path 'edpl_list.csv'

    """
    set -e


    gappa examine edpl \
    --jplace-path ${dedup_jplace} \
    --verbose \
    --threads ${task.cpus}
    """
}

process MakeSplit {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
        path sv_long_f
    output:
        path 'sv_multiplicity.csv'

"""
#!/usr/bin/env python
import csv

with open('${sv_long_f}', 'rt') as in_h, open('sv_multiplicity.csv', 'wt') as out_h:
    svl = csv.DictReader(in_h)
    t_svl = csv.writer(
        out_h,
    )
    for r in svl:
        t_svl.writerow([
            r['sv'],
            r['specimen'],
            r['count']
        ])
"""
}

process GappaSplit {
    container = "${container__gappa}"
    label = 'multithreaded'
    publishDir "${params.output}/placement/", mode: 'copy'

    input:
        path dedup_jplace
        path split_csv
    
    output:
        path 'specimen_jplace/*.jplace.gz'

    """
    set -e

    mkdir specimen_jplace

    gappa edit split \
    --jplace-path ${dedup_jplace} \
    --split-file ${split_csv} \
    --compress \
    --verbose \
    --threads ${task.cpus} \
    --out-dir specimen_jplace
    """
}

process Gappa_KRD {
    container = "${container__gappa}"
    label = 'mem_veryhigh'
    publishDir "${params.output}/placement/", mode: 'copy'
    errorStrategy 'finish'

    input:
        path specimen_jplace
    
    output:
        path 'krd/krd_matrix.csv.gz'

    """
    set -e

    gappa analyze krd \
    --jplace-path ${specimen_jplace} \
    --krd-out-dir krd/ \
    --krd-compress \
    --verbose \
    --threads ${task.cpus}

    """
}

process Gappa_ePCA {
    container = "${container__gappa}"
    label = 'mem_veryhigh'
    publishDir "${params.output}/placement/", mode: 'copy'
    errorStrategy 'finish'

    input:
        path specimen_jplace
    
    output:
        path 'ePCA/projection.csv'
        path 'ePCA/transformation.csv'

    """
    set -e

    gappa analyze edgepca \
    --jplace-path ${specimen_jplace} \
    --out-dir ePCA/ \
    --verbose \
    --threads ${task.cpus}

    ls -l ePCA

    """
}

process PplacerADCL {
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

process PplacerEDPL {
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

process PplacerPCA {
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
    tar xzvf ${refpkg_tgz_f} --no-overwrite-dir -C refpkg/ &&
    guppy epca ${dedup_jplace_f}:${sv_map_f} -c refpkg/ --out-dir pca/ --prefix epca &&
    guppy lpca ${dedup_jplace_f}:${sv_map_f} -c refpkg/ --out-dir pca/ --prefix lpca
    """
}

process PplacerAlphaDiversity {
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


process PplacerKR {
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
    tar xzvf ${refpkg_tgz_f} --no-overwrite-dir -C refpkg/
    guppy kr --list-out -c refpkg/ ${dedup_jplace_f}:${sv_map_f} |
    gzip > kr_distance.csv.gz
    """
}

process ClassifyDB_Prep {
    container = "${container__pplacer}"
    label = 'io_limited'
    afterScript "rm -r refpkg/"
    cache = false

    input:
        file refpkg_tgz_f
        file sv_map_f
    
    output:
        file 'classify.prep.db'
    

    """
    mkdir -p refpkg/
    tar xzvf ${refpkg_tgz_f} --no-overwrite-dir -C refpkg/
    rppr prep_db -c refpkg/ --sqlite classify.prep.db
    (echo "name,specimen"; cat ${sv_map_f}) |
    csvsql --table seq_info --insert --snifflimit 1000 --db sqlite:///classify.prep.db
    """
}

process ClassifySV {
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
        file 'classify.classified.db'

    """
    mkdir -p refpkg/
    tar xzvf ${refpkg_tgz_f} --no-overwrite-dir -C refpkg/
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

process ClassifyMCC {
    container = "${container__pplacer}"
    label = 'io_limited'
    cache = false
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        file classifyDB_classified
        file sv_weights_f

    output:
        file 'classify.mcc.db'

    """
    multiclass_concat.py -k \
    --dedup-info ${sv_weights_f} ${classifyDB_classified}
    cp ${classifyDB_classified} classify.mcc.db
    """
}

process ClassifyTables {
    container = "${container__pplacer}"
    label = 'io_limited'
    publishDir "${params.output}/classify", mode: 'copy'

    input:
        tuple val(rank), file(classifyDB_mcc), file(sv_map_for_tables_f)

    output:
        tuple val(rank), file("tables/by_specimen.${rank}.csv"), file("tables/by_taxon.${rank}.csv"), file("tables/tallies_wide.${rank}.csv")

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

process SharetableToMapWeight {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
        file (sharetable)

    output:
        file ("sv_sp_map.csv")
        file ("sv_weights.csv")

"""
#!/usr/bin/env python
import csv

sp_count = {}
with open('${sharetable}', 'rt') as st_h:
    st_r = csv.reader(st_h, delimiter='\\t')
    header = next(st_r)
    sv_name = header[3:]
    for r in st_r:
        sp_count[r[0]] = [int(c) for c in r[3:]]
weightsL = []
mapL = []
for sv_i, sv in enumerate(sv_name):
    sv_counts = [
        (sp, c[sv_i]) for sp, c in sp_count.items()
        if c[sv_i] > 0
    ]
    weightsL += [
        (sv, "{}__{}".format(sv, sp), c)
        for sp, c in 
        sv_counts
    ]
    mapL += [
        ("{}__{}".format(sv, sp), sp)
        for sp, c in 
        sv_counts
    ]
with open("sv_sp_map.csv", "w") as map_h:
    map_w = csv.writer(map_h)
    map_w.writerows(mapL)
with open("sv_weights.csv", "w") as weights_h:
    weights_w = csv.writer(weights_h)
    weights_w.writerows(weightsL) 
"""
}

process Dada2_convert_output {
    container "${container__dada2pplacer}"
    label 'io_mem'
    publishDir "${params.output}/sv/", mode: 'copy'
    errorStrategy "retry"

    input:
        file(final_seqtab_csv)

    output:
        file "dada2.sv.fasta"
        file "dada2.sv.map.csv"
        file "dada2.sv.weights.csv"

    """
    dada2-seqtab-to-pplacer \
    -s ${final_seqtab_csv} \
    -f dada2.sv.fasta \
    -m dada2.sv.map.csv \
    -w dada2.sv.weights.csv \
    """
}

process Extract_Taxonomy {
    container "${container__dada2pplacer}"
    label 'io_mem'
    publishDir "${params.output}/classify", mode: 'copy'
    errorStrategy "ignore"

    input:
        file (weights_csv)
        file (taxonomy_db)

    output:
        file "sv_taxonomy.csv"

"""
#!/usr/bin/env python
import csv
import pandas as pd
import sqlite3

sv = {
    r[0] for r in 
    csv.reader(open(
        "${weights_csv}",
        'rt'
    ))
}

tax_db_conn = sqlite3.connect("${taxonomy_db}")
tax_db_cur = tax_db_conn.cursor()
tax_db_cur.execute("CREATE INDEX IF NOT EXISTS idx_mcc_name_tax  ON multiclass_concat(name, tax_id)")

sv_classification = pd.concat([
    pd.read_sql("SELECT name, want_rank, taxa.tax_id, taxa.tax_name, taxa.rank, likelihood FROM multiclass_concat JOIN taxa ON multiclass_concat.tax_id = taxa.tax_id WHERE name=(?)",
        con = tax_db_conn,
        params=(
            sv,
        ))    
    for sv in sv
], ignore_index=True)

sv_classification.to_csv("sv_taxonomy.csv", index=False)

"""
}


//
//      Components for running the module independently
//

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/maliampi <ARGUMENTS>
    
    Required Arguments:
        --refpkg              Reference Package in tar.gz format
                        AND
    ONE of the following (for sv-counts-per-specimen):
        --sv_fasta            Fasta file with all the sequence variants to be placed
        --weights             Headerless CSV file with weights (shared_sv_id, specimen_sv_id, count)
        --map                 (specimen_sv_id, specimen)
                        OR
        --sv_fasta            Fasta file with all the sequence variants to be placed
        --sharetable          Mothur-style sharetable
                        OR
        --seqtable            DADA2 style seqtable
    Options:
      Common to all:
        --output              Directory to place outputs (default invocation dir)
                                Maliampi will create a directory structure under this directory
        -w                    Working directory. Defaults to `./work`
        -resume                 Attempt to restart from a prior run, only completely changed steps

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

workflow {
    if (
        params.help ||
        params.refpkg == null
    ) {
        helpMessage()
        exit 0
    }

    refpkg_tgz_f = file(params.refpkg)

    if (
        (params.sv_fasta != null) &&
        (params.weights != null ) &&
        (params.map != null)
    ) {
        map_f =  Channel.from( file(params.map) )
        weights_f = Channel.from( file(params.weights) )
        sv_fasta_f = Channel.from( file(params.sv_fasta) )
    }
    else if (
        (params.sv_fasta != null) &&
        (params.sharetable != null)
    ) {
        sv_fasta_f = file(params.sv_fasta)
        SharetableToMapWeight (
            file(params.sharetable)
        )
        map_f = SharetableToMapWeight.out[0]
        weights_f = SharetableToMapWeight.out[1]
        
    }
    else if (params.seqtable != null) {
        Dada2_convert_output(file(params.seqtable))
        sv_fasta_f = Dada2_convert_output.out[0]
        map_f = Dada2_convert_output.out[1]
        weights_f = Dada2_convert_output.out[2]
    } else {
        helpMessage()
        exit 0
    }

    place_classify_wf (
        sv_fasta_f,
        refpkg_tgz_f,
        weights_f,
        map_f,
        sv_long_f,
    )
}