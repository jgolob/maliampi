//
//  ASV via swarm 3. 
//    Useful particularly for newer illumina reads with faked binned qual scores :/
//
nextflow.enable.dsl=2

container__swarm = "quay.io/biocontainers/swarm:3.1.2--h9f5acd7_0"
container__vsearch = "quay.io/biocontainers/vsearch:2.22.1--hf1761c0_0"
container__fastatools = "golob/fastatools:0.8.5A"



// parameters for individual operation
// Defaults for parameters

// common
params.output = '.'
params.help = false

// dada2-sv
params.trimLeft = 0
params.maxN = 0
params.maxEE = 'Inf'
params.truncLenF = 0
params.truncLenR = 0
params.truncQ = 2
params.d = 1
params.chimera_method = 'consensus'


workflow swarm_wf {
    // These should be preprocessed through TrimGalore, etc
    take: miseq_pe_ch
    take: miseq_se_ch
    take: pyro_ch

    main:
    //
    // STEP 1: Merge Paired reads (by specimen)
    //
    MergePairs(
        miseq_pe_ch.map { [it[0], it[2], it[3]]}
    )

    MergePairs.out.mix(
        miseq_se_ch.map{ [it[0], it[2]]}
    ).mix(
        pyro_ch.map{[it[0], it[2]]}
    ).set{
        to_filter_ch
    }
    //
    // STEP 2: Quality filtering and trimming (by specimen)
    //
    FilterAndTrim(
        to_filter_ch
    )
    //
    // STEP 3: Dereplicate and relabel (by specimen). Sequence-based counts
    //
    SpecimenDereplicate(
        FilterAndTrim.out
    )
    //
    // STEP 4: Dereplicate across specimens and rename into ASVs
    //
    MakeDSVforSwarm(
        SpecimenDereplicate.out
            .map{
                it[1]
            }
            .toList()
    )

    //
    // STEP 5: Run Swarm
    //

    Swarm(
        MakeDSVforSwarm.out.dsv_fasta
    )

    //
    // STEP 6: Chimera / singleton removal
    //

    RemoveGlobalSingleton(
        Swarm.out.swarm_seeds_fasta
    )

    ChimeraRemoval(
        RemoveGlobalSingleton.out
    )

    //
    // Step 7: Convert to ASVs
    // 
    SwarmToASV(
        ChimeraRemoval.out,
        Swarm.out.swarm_clusters,
        MakeDSVforSwarm.out.sp_dsv_count
    )

    //
    //  STEP 8: Convert outputs
    //
    ConvertOutputs(
        SwarmToASV.out.asv_fasta,
        SwarmToASV.out.sp_asv_long
    )

    emit:
       sv_fasta         = SwarmToASV.out.asv_fasta
       sv_map           = ConvertOutputs.out.map
       sv_weights       = ConvertOutputs.out.weights
       sv_long          = SwarmToASV.out.sp_asv_long
       sv_sharetable    = ConvertOutputs.out.sharetable
       sv_table         = ConvertOutputs.out.seqtable
    // */
}

process MergePairs {
    container "${container__vsearch}"
    label 'io_limited'
    errorStrategy "ignore"

    input:
        tuple val(specimen), file(R1), file(R2)
    
    output:
        tuple val(specimen), file("${specimen}.merged.fastq.gz")
    """
    vsearch --fastq_mergepairs \
    ${R1} --reverse ${R2} \
    --fastqout "${specimen}.merged.fastq" \
    --threads ${task.cpus} \
    --fastq_eeout

    gzip "${specimen}.merged.fastq"
    """
}

process FilterAndTrim {
    container "${container__vsearch}"
    label 'io_limited'
    errorStrategy "ignore"

    input:
        tuple val(specimen), file(R1)
    
    output:
        tuple val(specimen), file("${specimen}.filtered.fasta.gz")
    """
    vsearch --fastq_filter \
    ${R1} \
    --fastq_maxee ${params.maxEE} \
    --fastq_maxns ${params.maxN} \
    --fastq_truncqual ${params.truncQ} \
    --fastq_stripleft ${params.truncLenF} \
    --fastq_stripright ${params.truncLenR} \
    --threads ${task.cpus} \
    --fastaout "${specimen}.filtered.fasta"

    gzip "${specimen}.filtered.fasta"
    """
}

process SpecimenDereplicate {
    container "${container__vsearch}"
    label 'multithread'
    errorStrategy "ignore"

    input:
        tuple val(specimen), file(R1)
    
    output:
        tuple val(specimen), file("${specimen}.derep.fasta.gz")
    """
    vsearch --derep_fulllength \
    ${R1} \
    --strand plus \
    --sizeout \
    --relabel ${specimen}. \
    --fasta_width 0 \
    --threads ${task.cpus} \
    --output "${specimen}.derep.fasta"

    gzip "${specimen}.derep.fasta"
    """
}

process MakeDSVforSwarm {
    container = "${container__fastatools}"
    label = 'io_mem'
    
    input:
        path specimen_fasta

    output:
        path 'dsv.fasta.gz', emit: dsv_fasta
        path 'sp_dsv_count.csv.gz', emit: sp_dsv_count

"""
#!/usr/bin/env python
import gzip
import fastalite
from collections import defaultdict
import re
import csv

re_id = re.compile(r'^(?P<read>.+);size=(?P<num>\\d+)')

fasta_fns = "${specimen_fasta}".split()
# A list to store specimen, sequence, count in long format
sp_seq_c = []

for fn in fasta_fns:
    specimen = fn.replace('.derep.fasta.gz', '')
    sp_seq_c += [
        (
            specimen,
            sr.seq,
            int(re_id.match(sr.id)['num']),
        )
        for sr in
        fastalite.fastalite(
            gzip.open(fn, 'rt')
        )
    ]
# Sequence-counts 
seq_count = defaultdict(int)
for (sp, seq, count) in sp_seq_c:
    seq_count[seq] += count
# Convert to a sorted list by count
seq_count_l = sorted([
    (seq, count)
    for (seq, count) in seq_count.items()
], key=lambda v: -1*v[1])

seq_dsv = {}
# Assign DSV IDs and output fasta
with gzip.open('dsv.fasta.gz', 'wt') as out_dsv:
    for seq_i, (seq, count) in enumerate(seq_count_l):
        dsv = f'DSV{seq_i + 1:05d}'
        out_dsv.write(f'>{dsv};size={count}\\n{seq}\\n')
        seq_dsv[seq] = dsv

# Finally output our DSV <-> specimen <-> count table
with gzip.open('sp_dsv_count.csv.gz', 'wt') as sdc_h:
    sdc_w = csv.writer(sdc_h)
    sdc_w.writerow(['specimen', 'dsv', 'count'])
    for sp, seq, count in sp_seq_c:
        sdc_w.writerow([
            sp,
            seq_dsv.get(seq),
            count
        ])
"""
}

process Swarm {
    container "${container__swarm}"
    label 'multithread'
    errorStrategy "ignore"

    input:
        path (R1)
    
    output:
        path 'swarm_seeds.fasta.gz', emit: swarm_seeds_fasta
        path 'swarm_clusters.txt', emit: swarm_clusters
    """
    zcat ${R1} | \
    swarm \
    -d ${params.d} \
    -f -z \
    --threads ${task.cpus} \
    -o swarm_clusters.txt \
    -w swarm_seeds.fasta

    gzip swarm_seeds.fasta
    """
}

process RemoveGlobalSingleton {
    container = "${container__fastatools}"
    label = 'io_mem'
    
    input:
        path swarm_seeds_fasta

    output:
        path "swarm_seeds.no_singleton.fasta.gz"

"""
#!/usr/bin/env python
import gzip
import fastalite
import re

re_id = re.compile(r'^(?P<specimen>.+);size=(?P<num>\\d+)')

with gzip.open("swarm_seeds.no_singleton.fasta.gz", 'wt') as out_h:
    for sr in fastalite.fastalite(gzip.open('${swarm_seeds_fasta}', 'rt')):
        count = int(re_id.search(sr.id)['num'])
        if count > 1:
            out_h.write(f'>{sr.id}\\n{sr.seq}\\n')
"""
}

process ChimeraRemoval {
    container "${container__vsearch}"
    label 'io_limited'
    errorStrategy "ignore"
    publishDir "${params.output}/sv/", mode: 'copy'

    input:
        file(R1)
    
    output:
        file("swarm.seeds_nochimera.fasta.gz")
    """
    vsearch --uchime3_denovo \
    ${R1} \
    --threads ${task.cpus} \
    --nonchimeras "swarm.seeds_nochimera.fasta"

    gzip "swarm.seeds_nochimera.fasta"
    """
}

process SwarmToASV {
    container = "${container__fastatools}"
    label = 'io_mem'
    publishDir "${params.output}/sv/", mode: 'copy'

    input:
        path swarm_seeds_fasta
        path swarm_clusters
        path specimen_dsv_count

    output:
        path 'swarm_ASV.fasta.gz', emit: asv_fasta
        path 'sp_asv_long.csv.gz', emit: sp_asv_long

"""
#!/usr/bin/env python
import gzip
import fastalite
import csv
import pandas as pd

# First get a cluster index for each DSV

DSV_clusterIdx = {}
for c_i, line in enumerate(open('${swarm_clusters}', 'rt')):
    DSV_clusterIdx.update({
        dsv.split(';')[0]: c_i
        for dsv in line.split()
    })

# Then try to get the ASV sequence for each cluster
cluster_seq = {}
for sr in fastalite.fastalite(
    gzip.open("${swarm_seeds_fasta}", 'rt')
):
    # What is the 'representitive' DSV for this ASV
    rep_dsv = sr.id.split(';')[0]
    # Which cluster is this DSV in?
    asv_clusterIdx = DSV_clusterIdx.get(rep_dsv, -1)
    if asv_clusterIdx == -1:
        print(sr.id)
        continue
    # Implicit else
    cluster_seq[asv_clusterIdx] = sr

# Output the ASVs in fasta format
# Cache the cluster <-> ASV_id
cluster_ASVid = {}
with gzip.open('swarm_ASV.fasta.gz', 'wt') as ASV_h:
    n = 1
    for cluster_idx, cluster_sr in cluster_seq.items():
        asv_id = f'ASV{n:05d}'
        cluster_ASVid[cluster_idx] = asv_id
        ASV_h.write(f'>{asv_id}\\n{cluster_sr.seq}\\n')
        n += 1

# Now the per-specimen long format


sp_dsv_l = pd.read_csv('${specimen_dsv_count}')
sp_dsv_l['clusterIdx'] = sp_dsv_l.dsv.apply(DSV_clusterIdx.get)
sp_dsv_l['sv'] = sp_dsv_l.clusterIdx.apply(cluster_ASVid.get)
sp_asv_l = sp_dsv_l.groupby(['specimen', 'sv']).sum(
    numeric_only=True
).reset_index()[['specimen', 'sv', 'count']]

sp_asv_l.to_csv('sp_asv_long.csv.gz', index=None)

reads_kept = sp_asv_l['count'].sum()
reads_filtered = sp_dsv_l['count'].sum() - reads_kept

print("Kept", reads_kept)
print("Filtered", reads_filtered)
print("Percentage Kept", (reads_kept) / (reads_kept + reads_filtered) * 100)
"""
}

process ConvertOutputs {
    container = "${container__fastatools}"
    label = 'io_mem'
    publishDir "${params.output}/sv/", mode: 'copy'
    
    input:
        path asv_fasta
        path sp_asv_long

    output:
        path 'swarm_sv.seq_table.csv', emit: seqtable
        path 'swarm_sv.share_table.txt', emit: sharetable
        path 'swarm_sv.map.csv', emit: map
        path 'swarm_sv.weights.csv', emit: weights

"""
#!/usr/bin/env python
import gzip
import pandas as pd
import fastalite

asv_seq = {
    sr.id: sr.seq
    for sr in fastalite.fastalite(
        gzip.open('${asv_fasta}', 'rt')
    )
}

svl = pd.read_csv('${sp_asv_long}')

svw = svl.pivot(
    index='specimen',
    columns='sv',
    values='count'
).fillna(0).astype(int)
# Sharetable

st = pd.DataFrame(
    index=svw.index
)
st['label'] = st.index
st['group'] = 'swarm'
st['numsvs'] = len(svw.columns)
for c in svw.columns:
    st[c] = svw[c]

st.to_csv(
    'swarm_sv.share_table.txt',
    sep='\\t',
    index=None
)
# Seqtable
svw.rename(asv_seq, axis=1).to_csv(
    'swarm_sv.seq_table.csv'
)

# Map and weight
# decide on representitive specimens for each SV
sv_repSp = {
    sv: sp
    for (sv, sp) in
    svl.sort_values('count', ascending=False).groupby('sv').first().specimen.items()
}
# add specimen-specific SV
svl['sv_sp'] = [
    sv if sv_repSp[sv] == sp else
    f'{sv}__{sp}'
    for (sp, sv) in 
    zip(
        svl.specimen,
        svl.sv
    )
]
# Map is sv_sp, sp
svl[['sv_sp', 'specimen']].drop_duplicates().to_csv(
    'swarm_sv.map.csv',
    index=None,
    header=None
)
# Weights is global sv, sp_sv, count
svl[['sv', 'sv_sp', 'count']].drop_duplicates().to_csv(
    'swarm_sv.weights.csv',
    index=None,
    header=None
)
"""
}

//
// standalone workflow for module
//

include { read_manifest } from './manifest'
include { output_failed } from './preprocess' params (
    output: params.output
)
include { preprocess_wf } from './preprocess'
params.manifest = null
// Function which prints help message text
def helpMessage() {
    log.info"""
    Workflow to make sequence variants via swarm3 from a manifest of reads

    Usage:

    nextflow run jgolob/maliampi/dada2.nf <ARGUMENTS>
    
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

    """.stripIndent()
}

workflow {
    if (params.manifest == null) {
        helpMessage()
        exit 0
    }

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

    swarm_wf(
        preprocess_wf.out.miseq_pe,
        preprocess_wf.out.miseq_se,
        preprocess_wf.out.pyro
    )

    // */

}