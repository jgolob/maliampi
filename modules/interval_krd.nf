//
//  Take a new set of specimen-jplace and 
//  generate new KRD between each specimen-jplace and a collection of prior specimen-jplace
//  Do not repeat pairwise for old-old or new-new. just old-new pairs
//
nextflow.enable.dsl=2

// Parameters
params.old_jplace = false
params.new_jplace = false
params.output = '.'
container__krd = "golob/gappa:0.5"

workflow {
    
    old_jplace_files = Channel.fromPath(params.old_jplace+"/*.jplace.gz")
    new_jplace_files = Channel.fromPath(params.new_jplace+"/*.jplace.gz")

    Gappa_KRD_1t1(
        new_jplace_files,
        old_jplace_files.toList(),
    )
    CombineKRDLong(
        Gappa_KRD_1t1.out.toList()
    )
}

process Gappa_KRD_1t1 {
    container = "${container__krd}"
    label = 'io_limited'

    input:
        path new_jplace
        path old_jplaces, stageAs: "that/"
    output:
        path "v${new_jplace.name.replace('.jplace.gz', '')}.krd_long.csv"

    
    """
    interval_krd.py \
    --new-jplace ${new_jplace} \
    --old-jplaces that/
    """
}

process Gappa_KRD_1tm {
    container = "${container__krd}"
    label = 'io_limited'

    input:
        tuple path(this_jplace), path(those_jplaces, stageAs: "those/")
    output:
        path "v${this_jplace.name.replace('.jplace.gz', '')}.krd_long.csv"

    
    """
    interval_krd.py \
    --new-jplace ${this_jplace} \
    --old-jplaces those/
    """
}

process CombineKRDLong {
    container = "${container__krd}"
    label = 'io_limited'
    publishDir "${params.output}/placement", mode: 'copy'

    input:
        path krd_longs
    output:
        path "krd_long.csv.gz"
    
"""
#!/usr/bin/env python3
import csv
import gzip

krd_long_fn = "${krd_longs}".split()

with gzip.open("krd_long.csv.gz", 'wt') as out_h:
    w = csv.DictWriter(out_h, fieldnames=['specimen_1', 'specimen_2', 'krd'])
    w.writeheader()
    for fn in krd_long_fn:
        with open(fn, 'rt') as in_h:
            r = csv.DictReader(in_h)
            w.writerows([
                row for row in r
            ])
"""
}