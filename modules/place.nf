nextflow.enable.dsl=2

/*
 * Engine adapters only.  Taxonomy and statistics intentionally live in their
 * own modules so that a placement engine cannot select a downstream stack.
 */

process AlignSV {
    container "${params.container__infernal}"
    label 'mem_veryhigh'

    input:
    path sv_fasta
    path cm

    output:
    path 'sv.aln.sto', emit: stockholm
    path 'sv.aln.scores', emit: scores

    script:
    """
    cmalign --cpu ${task.cpus} --noprob --dnaout --mxsize ${params.cmalign_mxsize} \\
      --sfile sv.aln.scores -o sv.aln.sto ${cm} ${sv_fasta}
    """
}

process CombineAln_SV_refpkg {
    container "${params.container__easel}"
    label 'mem_veryhigh'

    input:
    path sv_aln_sto
    path refpkg_aln_sto

    output:
    path 'combined.sto', emit: stockholm

    script:
    """
    esl-alimerge --dna -o combined.sto ${sv_aln_sto} ${refpkg_aln_sto}
    """
}

process ConvertAlnToFasta {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'
    errorStrategy 'retry'

    input:
    path combined_sto

    output:
    path 'combined.fasta', emit: fasta

    script:
    """
    python3 - <<'PY'
    from Bio import AlignIO

    with open('combined.fasta', 'wt') as out_h:
        AlignIO.write(AlignIO.read('${combined_sto}', 'stockholm'), out_h, 'fasta')
    PY
    """
}

process EPAngPlacement {
    container "${params.container__epang}"
    label 'mem_veryhigh'
    publishDir "${params.output}/placement", mode: 'copy'

    input:
    path refpkg_aln_fasta
    path combined_aln_fasta
    path model
    path ref_tree

    output:
    path 'dedup.jplace', emit: dedup_jplace

    script:
    """
    set -e
    epa-ng --split ${refpkg_aln_fasta} ${combined_aln_fasta}
    model=\$(cat ${model})
    epa-ng -t ${ref_tree} -s reference.fasta -q query.fasta -m \$model \\
      -T ${task.cpus} --baseball-heur
    mv epa_result.jplace dedup.jplace
    """
}

process PplacerPlacement {
    container "${params.container__pplacer}"
    label 'mem_veryhigh'
    afterScript 'rm -rf refpkg/ || true'
    publishDir "${params.output}/placement", mode: 'copy'

    input:
    path combined_aln_sto
    path refpkg_tgz

    output:
    path 'dedup.jplace', emit: dedup_jplace

    script:
    """
    mkdir -p refpkg/
    tar xzvf ${refpkg_tgz} --no-overwrite-dir -C ./refpkg
    pplacer -p -j ${task.cpus} --inform-prior \\
      --prior-lower ${params.pplacer_prior_lower} --map-identity \\
      -c refpkg/ ${combined_aln_sto} -o dedup.jplace
    """
}

process GappaSplit {
    container "${params.container__gappa}"
    label 'multithread'
    publishDir "${params.output}/placement", mode: 'copy'

    input:
    path dedup_jplace
    path split_csv

    output:
    path 'specimen_jplace/**/*.jplace.gz', emit: specimen_jplaces

    script:
    """
    set -e
    mkdir specimen_jplace
    gappa edit split --jplace-path ${dedup_jplace} --split-file ${split_csv} \\
      --compress --verbose --threads ${task.cpus} --out-dir specimen_jplace
    """
}

process PlacementIdentity {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'
    publishDir "${params.output}/placement", mode: 'copy'

    input:
    path dedup_jplace
    path refpkg_identity
    path sv_h5ad
    val placer

    output:
    path 'placement.identity.json', emit: identity

    script:
    """
    python3 - <<'PY'
    import hashlib
    import json

    with open('${refpkg_identity}') as in_h:
        refpkg = json.load(in_h)
    with open('${dedup_jplace}') as in_h:
        jplace = json.load(in_h)
    tree = jplace.get('tree')
    if not isinstance(tree, str) or not tree:
        raise SystemExit('jplace does not contain a tree')
    payload = {
        'refpkg_id': refpkg['refpkg_id'],
        'tree_sha256': hashlib.sha256(tree.encode()).hexdigest(),
        'sv_artifact_sha256': hashlib.sha256(open('${sv_h5ad}', 'rb').read()).hexdigest(),
        'placer': '${placer}',
    }
    with open('placement.identity.json', 'w') as out_h:
        json.dump(payload, out_h, indent=2, sort_keys=True)
        out_h.write('\\n')
    PY
    """
}
