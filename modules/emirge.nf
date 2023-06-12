#!/usr/bin/env nextflow

// containers
container__emirge = "golob/emirge:0.62.1G"
container__fastcombineseqtab = "golob/dada2-fast-combineseqtab:0.5.0__1.12.0__BCW_0.3.1"
container__dada2pplacer = "golob/dada2-pplacer:0.8.0__bcw_0.3.1A"
container__trimgalore = 'quay.io/biocontainers/trim-galore:0.6.6--0'
container__fastatools = "golob/fastatools:0.8.0A"
// Using DSL-2
nextflow.enable.dsl=2


// Set default values for parameters
params.manifest = false
params.ref = false
params.output = '.'
params.help = false
params.nopreprocess = false

params.max_iterations = 10



// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run emirge.nf <ARGUMENTS>

    NOTE:   This script expects paired-end FASTQ data, and will not download any other type

    Required Arguments:
      --manifest            CSV format. Must have at least the following columns:
                                - specimen -> unique descriptor
                                - R1 -> path to first of paired end reads
                            Optionally will have
                                - R2 -> path to second of paired end reads
                                    (Can be empty)

      --ref                 FASTA file of 16s rRNA gene alleles to use as our reference.
                            Must be gzipped.

      --output              Folder to write output files, which will be organized to:
                                args.output/SRRxxxxxxxx.fasta.gz
                                Where the fasta file is the recovered 16S rRNA alleles
      
    Optional arguments:
      --nopreprocess        Skip trim step (only if data was previously preprocessed)

    Output Files:
    """.stripIndent()
}



workflow {
    main:
        // Show help message if the user specifies the --help flag at runtime
        if (params.help || params.manifest == false || params.output == false || params.ref == false){
            // Invoke the function above which prints the help message
            helpMessage()
            // Exit out and do not run anything else
            exit 0
        }

        // Make sure that --output ends with a trailing "/" character
        if (!params.output.endsWith("/")){
            output_folder = params.output.concat("/")
        } else {
            output_folder = params.output
        }
        
        manifest = Channel.fromPath(
            file(params.manifest)
        )
        .splitCsv(header: true)
        // Great. Split now into paired, single end, and invalid
        .branch {
            pe: (it.specimen != null) && (it.R1 != null) && (it.R1.trim() != '') && (it.R2 != null) && (it.R2.trim() != '')
            se: (it.specimen != null) && (it.R1 != null) && (it.R1.trim() != '')
            invalid: true
        }
    
        manifest.invalid.view()
        // Next check that the files are not empty and map
        pe_ch = manifest.pe.filter{
            (!file(it.R1).isEmpty()) && (!file(it.R2).isEmpty())
        }.map{[
            it.specimen.trim(),
            file(it.R1.trim()),
            file(it.R2.trim()),
        ]}

        se_ch = manifest.se.filter{
            (!file(it.R1).isEmpty())
        }.map{[
            it.specimen.trim(),
            file(it.R1.trim()),
        ]}

        // Trimgalore to clean up the reads...
        if (params.nopreprocess) {
            pe_ft_ch = pe_ch
            se_ft_ch = se_ch
        } else {
            TrimGalore(pe_ch).set{
                pe_ft_ch
            }
            TrimGaloreSE(se_ch).set{
                se_ft_ch
            }

        }
        
        


        // Get the max read length 
        MaxReadLen_SE(se_ft_ch)
        MaxReadLen_PE(pe_ft_ch)
  
        // For PE only, merge the pairs to get the inset lengths
        MergePairs(
            MaxReadLen_PE.out
        )
        // OK, if the read pairs FAIL to merge, go ahead and try to salvage as single end reads.
        pe_post_merge = MergePairs.out.branch{
            failed: it[4].strip() == '-nan' || it[5].strip() == '-nan'
            merged: true
        }
        
        // Convert the reference (fasta.gzip) into an index for use
        EMIRGE_MakeRef(
            file(params.ref)
        )
        
        // Then use the results to run emirge on these specimens
        EMIRGE_PE(
            pe_post_merge.merged,
            EMIRGE_MakeRef.out.ref_fasta,
            EMIRGE_MakeRef.out.ref_btindex
        )
      // Run EMIRGE on SE reads
        EMIRGE_SE(
            MaxReadLen_SE.out.mix(
                pe_post_merge.failed.map{
                    [it[0], it[1], it[3].strip()]
                }
            ),
            EMIRGE_MakeRef.out.ref_fasta,
            EMIRGE_MakeRef.out.ref_btindex
        )

    // Make seq_longs

        MakeSeqLong(
            EMIRGE_PE.out.mix(EMIRGE_SE.out).map{[
                it[0],
                it[1]
            ]}
        )
        
        CombineSeqLong(
            MakeSeqLong.out
                .map{ it[2]  }
		        .toList()
        ) 
// */
}

process MaxReadLen_SE {
    container = "${container__emirge}"
    label = "io_limited"
    //errorStrategy 'finish'

    input:
        tuple val(SRR), file(R1)
    
    output:
        tuple val(SRR), file(R1), stdout
    
    """
    #!/usr/bin/env python

    from Bio import SeqIO
    import gzip

    max_read_len = 0
    for sr in SeqIO.parse(gzip.open('${R1}', 'rt'), 'fastq'):
        max_read_len = max([
            max_read_len,
            len(sr)
        ])

    print(max_read_len)
    """
}

process MaxReadLen_PE {
    container = "${container__emirge}"
    label = "io_limited"
    //errorStrategy 'finish'

    input:
        tuple val(SRR), file(R1), file(R2)
    
    output:
        tuple val(SRR), file(R1), file(R2), stdout
    
    """
    #!/usr/bin/env python

    from Bio import SeqIO
    import gzip

    max_read_len = 0
    for sr in SeqIO.parse(gzip.open('${R1}', 'rt'), 'fastq'):
        max_read_len = max([
            max_read_len,
            len(sr)
        ])

    for sr in SeqIO.parse(gzip.open('${R2}', 'rt'), 'fastq'):
        max_read_len = max([
            max_read_len,
            len(sr)
        ])

    print(max_read_len)
    """
}

process MergePairs {
    container = "${container__emirge}"
    label = "multithread"

    input: 
        tuple val(SRR), file(R1), file(R2), val(MAX_LEN)
    output: 
        tuple val(SRR), file(R1), file(R2), val(MAX_LEN), env(INSERT_MEAN), env(INSERT_STD)
    
    """
    vsearch --fastq_mergepairs \
    ${R1} \
    --reverse ${R2} \
    --threads ${task.cpus} \
    --fastq_minovlen 1 \
    --fastaout /dev/null \
    >> vsearch.out.txt 2>&1
    
    cat vsearch.out.txt
    
    INSERT_MEAN=\$(cat vsearch.out.txt | grep 'Mean fragment length' | sed 's/  Mean.*//' | sed 's/\\..*//')
    INSERT_STD=\$(cat vsearch.out.txt | grep 'Standard deviation of fragment length' | sed 's/  Standard.*//' | sed 's/\\..*//')
    """
}

process EMIRGE_MakeRef {
    container = "${container__emirge}"
    label = "multithread"
    //publishDir "${params.output}/EMIRGE/${SRR}", mode: 'copy'
    maxForks 5
    errorStrategy 'finish'

    input:
        path(ref_fasta_gz)
    
    output:
        path('SSU_candidates.fasta'), emit: ref_fasta
        path('SSU_candidates_btindex.*'), emit: ref_btindex
    
    """
    gunzip -c ${ref_fasta_gz} > SSU_candidates.fasta
    bowtie-build SSU_candidates.fasta SSU_candidates_btindex --threads ${task.cpus}
    """
}

process EMIRGE_SE {
    container = "${container__emirge}"
    label = "multithread"
    publishDir "${params.output}/EMIRGE/${SRR}", mode: 'copy'
    errorStrategy 'ignore'
    maxForks 5

    input:
        tuple val(SRR), file(R1), val(MAX_READ_LEN)
        path (ref_fasta)
        path (ref_btindex)
    
    output:
        tuple val(SRR), file("${SRR}_final.fasta"), file("${SRR}_mapped_reads.sam"), file("${SRR}_bowtie.final.PE.bam"), file("${SRR}.emirge.log"), file("${SRR}_emirge.tgz")
    
    """
    emirge.py \
    ${SRR}_emirge_run \
    --processors ${task.cpus} \
    -1 ${R1} \
    -l ${MAX_READ_LEN.trim()} \
    -f ${ref_fasta} \
    -b ${ref_btindex[0].getSimpleName()}  \
    --phred33 \
    --iterations ${params.max_iterations} \
     >> ${SRR}.emirge.log 2>&1 || true 
    
    cd ${SRR}_emirge_run
    LAST_ITER=\$(ls -t iter.* | head -n 1 | sed 's/://')
        emirge_rename_fasta.py --no_N \$LAST_ITER > ../${SRR}_final.fasta
    ln -s \$LAST_ITER iter.final
    cp iter.final/bowtie.iter.??.PE.bam ../${SRR}_bowtie.final.PE.bam

    cd ..
    samtools view -F 4 ${SRR}_bowtie.final.PE.bam > ${SRR}_mapped_reads.sam
    tar czvf ${SRR}_emirge.tgz ${SRR}_emirge_run/

    rm -rf ${SRR}_emirge_run/iter.*
    """
}

process EMIRGE_PE {
    container = "${container__emirge}"
    label = "multithread"
    publishDir "${params.output}/EMIRGE/${SRR}", mode: 'copy'
    maxForks 5
    errorStrategy 'ignore'

    input:
        tuple val(SRR), file(R1), file(R2), val(MAX_READ_LEN), val(INSERT_MEAN), val(INSERT_STD)
        path (ref_fasta)
        path (ref_btindex)
    
    output:
        tuple val(SRR), file("${SRR}_final.fasta"), file("${SRR}_mapped_reads.sam"), file("${SRR}_bowtie.final.PE.bam"), file("${SRR}.emirge.log"), file("${SRR}_emirge.tgz")
    
    """
    gunzip -c ${R2} > ${SRR}.R2.fasta

    emirge.py \
    ${SRR}_emirge_run \
    --processors ${task.cpus} \
    -1 ${R1} \
    -l ${MAX_READ_LEN.trim()} \
    -2 ${SRR}.R2.fasta \
    -i ${INSERT_MEAN.trim()} \
    -s ${INSERT_STD.trim()} \
    -f ${ref_fasta} \
    -b ${ref_btindex[0].getSimpleName()} \
    --phred33 \
    --iterations ${params.max_iterations} \
     >> ${SRR}.emirge.log 2>&1 || true
    
    cd ${SRR}_emirge_run
    LAST_ITER=\$(ls -t iter.* | head -n 1 | sed 's/://')
    emirge_rename_fasta.py --no_N \$LAST_ITER > ../${SRR}_final.fasta
    ln -s \$LAST_ITER iter.final
    cp iter.final/bowtie.iter.??.PE.bam ../${SRR}_bowtie.final.PE.bam
    

    cd ..
    samtools view -F 4 ${SRR}_bowtie.final.PE.bam > ${SRR}_mapped_reads.sam
    tar czvf ${SRR}_emirge.tgz ${SRR}_emirge_run/
    rm -rf ${SRR}_emirge_run/iter.*
    """
}

process MakeSeqLong {
    container = "${container__fastatools}"
    label = "io_limited"
    errorStrategy 'ignore'

    input:
        tuple val(SRR), path(SRR_fasta)

    output:
        tuple val(SRR), path(SRR_fasta), path("${SRR}.seq_long.csv")

"""
#!/usr/bin/env python
import fastalite
import re
import csv

re_desc = re.compile(r'^(?P<id>^\\d+\\|\\w+) Prior=(?P<prior>\\d\\.\\d+) Length=(?P<length>\\d+) NormPrior=(?P<norm_prior>\\d\\.\\d+)\$')

seq_long = []

for sr in fastalite.fastalite(open('${SRR_fasta}', 'rt')):
    m = re_desc.search(sr.description)
    if m is None:
        print(sr.description)
        continue
    # Implicit else
    seq_long.append({
        'id': m['id'],
        'seq': sr.seq,
        'length': int(m['length']),
        'n_amb': len([b for b in sr.seq if b.upper() not in {'A', 'C', 'T', 'G'} ]),
        'fract': float(m['prior']),
        'fract_n': float(m['norm_prior']),
        'specimen': '${SRR}',
        
    })

min_fract = min([v['fract'] for v in seq_long])
if min_fract != 0:
	for v in seq_long:
	    v['count'] = int(v['fract'] / min_fract)
else:
	for v in seq_long:
		v['count'] = int(v['fract']*1000)

with open('${SRR}.seq_long.csv', 'wt') as out_h:
    w = csv.DictWriter(
        out_h,
        fieldnames=[
            'specimen',
            'seq',
            'count',
            'fract',
            'fract_n',
            'length',
            'n_amb',
            'id'
        ]
    )
    w.writeheader()
    w.writerows(seq_long)
"""
}

process CombineSeqLong {
    container = "${container__dada2pplacer}"
    label = "io_mem"
    publishDir "${params.output}/", mode: 'copy'

    input:
        file(seq_longs)
    output:
        path("combined.seq_long.csv")
        path("combined.sv.fasta")

"""
#!/usr/bin/env python
import pandas as pd
csl = pd.concat([
    pd.read_csv(fn)
    for fn in "${seq_longs}".split()
], ignore_index=True)
seq_sv = {
    seq: 'sv_{:05d}'.format(s_i+1)
    for s_i, seq in 
    enumerate(
        csl.groupby('seq').sum()['count'].sort_values(ascending=False).index
    )
}
csl['sv'] = csl.seq.apply(seq_sv.get)
csl.to_csv("combined.seq_long.csv", index=None)

with open('combined.sv.fasta', 'wt') as sv_out:
    for seq, sv in seq_sv.items():
        sv_out.write(">{}\\n{}\\n".format(
            sv,
            seq
        ))
"""
}

process TrimGaloreSE {
    container "${container__trimgalore}"
    label 'io_limited'
    errorStrategy 'ignore'

    input:
    tuple val(specimen), file(R1)

    output:
    tuple val(specimen), file("${specimen}.R1.tg.fastq.gz")

    """
    set -e

    cp ${R1} R1.fastq.gz

    trim_galore \
    --gzip \
    --cores ${task.cpus} \
    R1.fastq.gz

    rm R1.fastq.gz
    mv R1_trimmed.fq.gz "${specimen}.R1.tg.fastq.gz"
    """
}

process TrimGalore {
    container "${container__trimgalore}"
    label 'io_limited'
    errorStrategy 'ignore'

    input:
    tuple val(specimen), file(R1), file(R2)

    output:
    tuple val(specimen), file("${specimen}.R1.tg.fastq.gz"), file("${specimen}.R2.tg.fastq.gz")

    """
    set -e

    cp ${R1} R1.fastq.gz
    cp ${R2} R2.fastq.gz

    trim_galore \
    --gzip \
    --cores ${task.cpus} \
    --paired \
    R1.fastq.gz R2.fastq.gz

    rm R1.fastq.gz
    rm R2.fastq.gz
    mv R1_val_1.fq.gz "${specimen}.R1.tg.fastq.gz"
    mv R2_val_2.fq.gz "${specimen}.R2.tg.fastq.gz"
    """
}