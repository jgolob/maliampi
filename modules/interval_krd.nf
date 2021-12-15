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
container__gappa = "golob/gappa:0.2"

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
    container = "${container__gappa}"
    label = 'multithread'

    input:
        path new_jplace
        path old_jplaces
    output:
        path "v${new_jplace.name.replace('.jplace.gz', '')}.krd_long.csv"

    
    """
    interval_krd.py \
    --new-jplace ${new_jplace} \
    --old-jplaces "${old_jplaces}"
    """
}

process CombineKRDLong {
    container = "${container__gappa}"
    label = 'io_limited'
    publishDir "${params.output}/", mode: 'copy'

    input:
        path krd_longs
    output:
        path "interval_krd_long.csv.gz"
    
"""
#!/usr/bin/env python3
import csv
import gzip

krd_long_fn = "${krd_longs}".split()

with gzip.open("interval_krd_long.csv.gz", 'wt') as out_h:
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