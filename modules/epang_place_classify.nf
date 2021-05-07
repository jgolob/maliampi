//
//  EPA-ng Place and Classify
//
nextflow.enable.dsl=2

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
container__infernal = "quay.io/biocontainers/infernal:1.1.4--h779adbc_0"
container__fastatools = "golob/fastatools:0.8.0A"
container__pplacer = "golob/pplacer:1.1alpha19rc_BCW_0.3.1A"
container__dada2pplacer = "golob/dada2-pplacer:0.8.0__bcw_0.3.1A"
container__easel = 'quay.io/biocontainers/easel:0.47--h516909a_0'
container__epang = "quay.io/biocontainers/epa-ng:0.3.8--h9a82719_1"
container__gappa = 'quay.io/biocontainers/gappa:0.7.1--h9a82719_1'


// includes
include { Dada2_convert_output } from './dada2' params (
    output: params.output
)

workflow epang_place_classify_wf {
    take:
        sv_fasta_f
        refpkg_tgz_f
        sv_weights_f
        sv_map_f
        sv_long_f

    main:

    //
    //  PLACEMENT
    //


    // Step 0. Extract bits from the reference
    ExtractRefpkg(
        refpkg_tgz_f
    )

    //
    // Step 1. Align the SV
    // 
    AlignSV(
        sv_fasta_f,
        ExtractRefpkg.out.cm
    )

    //
    //  Step 2. Combine SV and refpkg alignments
    //
    CombineAln_SV_refpkg(
        AlignSV.out[0],
        ExtractRefpkg.out.ref_aln_sto
    )

    // Step 2c. Convert combined alignment to FASTA format
    ConvertAlnToFasta(
        CombineAln_SV_refpkg.out
    )

    //
    //  Step 3. Place SV via epa-ng
    //

    EPAngPlacement(
        ExtractRefpkg.out.ref_aln_fasta,
        ConvertAlnToFasta.out,
        ExtractRefpkg.out.model,
        ExtractRefpkg.out.tree
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
        ExtractRefpkg.out.leaf_info,
        ExtractRefpkg.out.taxonomy,
    )

    Gappa_Classify(
        EPAngPlacement.out,
        MakeEPAngTaxonomy.out
    )

    Gappa_Extract_Taxonomy(
        Gappa_Classify.out,
        ExtractRefpkg.out.taxonomy
    )
    want_ranks = Channel.from(
        'species',
        'genus',
        'family',
        'class',
        'order',
        'phylum'
    )
    Make_Wide_Tax_Table(
        sv_long_f,
        Gappa_Extract_Taxonomy.out[0],
        want_ranks
    )

    emit:
        jplace_dedup = EPAngPlacement.out
        taxonomy = Gappa_Extract_Taxonomy.out[0]

}
process AlignSV {
    container = "${container__infernal}"
    label = 'mem_veryhigh'

    input:
        path sv_fasta_f
        path cm
    
    output:
        path "sv.aln.sto"
        path "sv.aln.scores"
        
    
    """
    cmalign \
    --cpu ${task.cpus} --noprob --dnaout --mxsize ${params.cmalign_mxsize} \
    --sfile sv.aln.scores -o sv.aln.sto \
    ${cm} ${sv_fasta_f}
    """
}




process CombineAln_SV_refpkg {
    container = "${container__easel}"
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

process ExtractRefpkg {
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
        path 'refpkg.cm', emit: cm

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

with open('refpkg.cm', 'wb') as cm_h:
    cm_h.write(
        tar_h.extractfile(
            tar_contents_dict[contents['files'].get('profile')]
        ).read()
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
if 'raxml_ng_model' in contents['files']:
    with open('model.txt', 'wt') as out_h:
        out_h.write(tar_h.extractfile(
            tar_contents_dict[contents['files'].get('raxml_ng_model')]
        ).read().decode('utf-8'))

else:
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
    errorStrategy 'ignore'

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

process Make_Wide_Tax_Table {
    container "${container__dada2pplacer}"
    label 'io_mem'
    publishDir "${params.output}/classify", mode: 'copy'
    //errorStrategy "ignore"

    input:
        path sv_long
        path sv_taxonomy
        val want_rank

    output:
        path "tables/taxon_wide_ra.${want_rank}.csv", emit: ra
        path "tables/taxon_wide_nreads.${want_rank}.csv", emit: nreads

"""
#!/usr/bin/env python
import pandas as pd
import os

try:
    os.makedirs('tables')
except:
    pass

sv_long = pd.read_csv("${sv_long}").rename({
    'count': 'nreads'
}, axis=1)
# Add in rel abund
for sp, sp_sv in sv_long.groupby('specimen'):
    sv_long.loc[sp_sv.index, 'fract'] = sp_sv.nreads / sp_sv.nreads.sum()

sv_taxonomy = pd.read_csv('${sv_taxonomy}')
sv_long_tax = pd.merge(
    sv_long,
    sv_taxonomy[sv_taxonomy.want_rank == '${want_rank}'],
    on='sv',
    how='left'
)

sp_tax = sv_long_tax.groupby(['specimen', 'tax_name']).sum().reset_index()[[
    'specimen',
    'tax_name',
    'nreads',
    'fract'
]]

sp_tax_wide_ra = sp_tax.pivot(
    index='specimen',
    columns='tax_name',
    values='fract'
).fillna(0)
# Sort by mean RA
sp_tax_wide_ra = sp_tax_wide_ra[sp_tax_wide_ra.mean().sort_values(ascending=False).index]

sp_tax_wide_ra.to_csv("tables/taxon_wide_ra.${want_rank}.csv")

sp_tax_wide_nreads = sp_tax.pivot(
    index='specimen',
    columns='tax_name',
    values='nreads'
).fillna(0)
# Sort by mean RA
sp_tax_wide_nreads = sp_tax_wide_nreads[sp_tax_wide_ra.mean().sort_values(ascending=False).index].astype(int)
sp_tax_wide_nreads.to_csv("tables/taxon_wide_nreads.${want_rank}.csv")

"""
}

process Gappa_Extract_Taxonomy {
    container "${container__dada2pplacer}"
    label 'io_mem'
    publishDir "${params.output}/classify", mode: 'copy'
    //errorStrategy "ignore"

    input:
        path gappa_taxonomy
        path refpkg_taxtable

    output:
        path "sv_taxonomy.csv"
        path refpkg_taxtable

"""
#!/usr/bin/env python
import pandas as pd

MIN_AFRACT = 0
RANKS = [
    'superkingdom',
    'phylum',
    'class',
    'order',
    'family',
    'genus',
    'species',
]
RANK_DEPTH = {
    i+1: r for (i, r) in enumerate(RANKS)
}
refpkg_taxtable = pd.read_csv("${refpkg_taxtable}")
tax_name_to_id = {
    row.tax_name: row.tax_id for
    idx, row in refpkg_taxtable.iterrows()
}
epa_tax = pd.read_csv('${gappa_taxonomy}', sep='\t')
epa_tax['lineage']=epa_tax.taxopath.apply(lambda tp: tp.split(';'))
epa_tax['rank_depth']=epa_tax.lineage.apply(len)
sv_tax_list = []
for sv, sv_c in epa_tax[epa_tax.taxopath != 'DISTANT'].groupby('name'):
    sv_tax = pd.DataFrame()
    rank = None
    tax_name = None
    lineage = None
    afract = None    
    for rank_depth, want_rank in RANK_DEPTH.items():
        sv_depth = sv_c[sv_c.rank_depth == rank_depth]
        if len(sv_depth) > 0 and sv_depth.afract.sum() >= MIN_AFRACT:
            # Something at this depth, and the cumulative fract likelihood is above our threshold
            rank = want_rank
            tax_name = " / ".join(sv_depth.lineage.apply(lambda L: L[-1]))
            ncbi_tax_id = ",".join([str(tax_name_to_id.get(tn, -1)) for tn in sv_depth.lineage.apply(lambda L: L[-1])])
            lineage = ";".join(sv_depth.lineage.iloc[0][:-1] + [tax_name])
            afract = sv_depth.afract.sum()
            
            
        sv_tax.loc[rank, 'sv'] = sv
        sv_tax.loc[rank, 'want_rank'] = want_rank
        sv_tax.loc[rank, 'rank'] = rank
        sv_tax.loc[rank, 'rank_depth'] = rank_depth
        sv_tax.loc[rank, 'tax_name'] = tax_name
        sv_tax.loc[rank, 'ncbi_tax_id'] = ncbi_tax_id
        sv_tax.loc[rank, 'lineage'] = lineage
        sv_tax.loc[rank, 'afract'] = afract
        sv_tax.loc[rank, 'ambiguous'] = len(sv_depth) != 1
    sv_tax_list.append(sv_tax)


sv_taxonomy = pd.concat(sv_tax_list, ignore_index=True)
sv_taxonomy['rank_depth'] = sv_taxonomy.rank_depth.astype(int)
sv_taxonomy.to_csv('sv_taxonomy.csv', index=None)

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
    errorStrategy 'ignore'

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
    errorStrategy 'ignore'

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

process WeightMaptoLong {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
        path (weight)
        path (map)

    output:
        path ("sp_sv_long.csv")

"""
#!/usr/bin/env python
import csv

specimen_comSV = {
    r[0]: r[1]
    for r in 
    csv.reader(
        open('${map}', 'rt')
    )
}
with open('${weight}', 'rt') as w_h, open("sp_sv_long.csv", 'wt') as sv_long_h:
    w_r = csv.reader(w_h)
    svl_w = csv.writer(sv_long_h)
    svl_w.writerow((
        'specimen',
        'sv',
        'count'
    ))    
    for row in w_r:
        svl_w.writerow((
            specimen_comSV[row[1]],
            row[0],
            int(row[2])
        ))
    
"""
}

process SharetableToMapWeight {
    container = "${container__fastatools}"
    label = 'io_limited'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
        path (sharetable)

    output:
        path ("sv_sp_map.csv"), emit: sv_map
        path ("sv_weights.csv"), emit: sv_weights
        path ("sp_sv_long.csv"), emit: sp_sv_long

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
sv_long = []
for sv_i, sv in enumerate(sv_name):
    sv_counts = [
        (sp, c[sv_i]) for sp, c in sp_count.items()
        if c[sv_i] > 0
    ]
    if len(sv_counts) == 0:
        continue
    # Implicit else
    shared_sv = "{}:{}".format(sv, sorted(sv_counts, key=lambda v: -1*v[1])[0][0])
    sv_long += [
        (sp, shared_sv, c[sv_i]) for sp, c in sp_count.items()
        if c[sv_i] > 0
    ]    
    weightsL += [
        (shared_sv, "{}:{}".format(sv, sp), c)
        for sp, c in 
        sv_counts
    ]
    mapL += [
        ("{}:{}".format(sv, sp), sp)
        for sp, c in 
        sv_counts
    ]
with open("sv_sp_map.csv", "w") as map_h:
    map_w = csv.writer(map_h)
    map_w.writerows(mapL)
with open("sv_weights.csv", "w") as weights_h:
    weights_w = csv.writer(weights_h)
    weights_w.writerows(weightsL)
with open("sp_sv_long.csv", 'wt') as svl_h:
    svl_w = csv.writer(svl_h)
    svl_w.writerow((
        'specimen',
        'sv',
        'count'
    ))
    svl_w.writerows(sv_long)
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
params.sv_fasta = null
params.weights = null
params.map = null
params.sharetable = null
params.seqtable = null
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
        WeightMaptoLong(
            weights_f,
            map_f,
        )
        sv_long_f = WeightMaptoLong.out
    }
    else if (
        (params.sv_fasta != null) &&
        (params.sharetable != null)
    ) {
        sv_fasta_f = file(params.sv_fasta)
        SharetableToMapWeight (
            file(params.sharetable)
        )
        map_f = SharetableToMapWeight.out.sv_map
        weights_f = SharetableToMapWeight.out.sv_weights
        sv_long_f = SharetableToMapWeight.out.sp_sv_long
        
    }
    else if (params.seqtable != null) {
        Dada2_convert_output(file(params.seqtable))
        sv_fasta_f = Dada2_convert_output.out[0]
        map_f = Dada2_convert_output.out[1]
        weights_f = Dada2_convert_output.out[2]
        sv_long_f = Dada2_convert_output.out.sv_long
    } else {
        helpMessage()
        exit 0
    }

    epang_place_classify_wf (
        sv_fasta_f,
        refpkg_tgz_f,
        weights_f,
        map_f,
        sv_long_f,
    )
}