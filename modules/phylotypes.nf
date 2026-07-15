nextflow.enable.dsl=2

process Phylotypes {
    container "${params.container__phylotypes}"
    // The independent image intentionally defaults to the phylotypes CLI.
    // Clear it so Nextflow's own `/bin/bash -c …` task command runs normally.
    containerOptions "--entrypoint ''"
    label 'mem_veryhigh'

    input:
    tuple val(threshold), val(threshold_token), path(dedup_jplace)
    val distance
    val lwr_overlap

    output:
    tuple val(threshold), val(threshold_token), path("threshold-${threshold_token}/sv_phylotype.csv"), emit: external_mapping

    script:
    """
    mkdir threshold-${threshold_token}
    phylotypes --jplace ${dedup_jplace} --out threshold-${threshold_token}/sv_phylotype.csv \\
      --threshold_pd ${threshold} --distance ${distance} --lwr-overlap ${lwr_overlap}
    """
}

process ImportPhylotypeMapping {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_mem'
    publishDir "${params.output}/phylotypes", mode: 'copy'

    input:
    tuple val(threshold), val(threshold_token), path(mapping)

    output:
    tuple val(threshold), val(threshold_token), path("phylotype.${threshold_token}.parquet"), emit: mapping

    script:
    """
    maliampi-phylotype-mapping ${mapping} phylotype.${threshold_token}.parquet
    """
}

process AggregatePhylotype {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_mem'
    publishDir "${params.output}/phylotypes", mode: 'copy'

    input:
    tuple val(threshold), val(threshold_token), path(mapping), path(sv_h5ad), path(placement_identity)

    output:
    tuple val(threshold), val(threshold_token), path("phylotype.${threshold_token}.h5ad"), emit: abundance_h5ad

    script:
    """
    maliampi-aggregate ${sv_h5ad} ${mapping} phylotype.${threshold_token}.h5ad \
      --identity-column phylotype --metadata-json '{"threshold": ${threshold}}' \\
      --placement-identity ${placement_identity}
    """
}

process BuildPhylotypeSet {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'
    publishDir "${params.output}/phylotypes", mode: 'copy'

    input:
    path dedup_jplace
    path placement_identity
    path mappings
    val dataset_id
    val thresholds
    val distance
    val lwr_overlap

    output:
    path 'phylotype-set.tar.gz', emit: set_bundle

    script:
    """
    python3 - <<'PY'
    import hashlib
    import json
    import pathlib
    import shutil
    import tarfile

    root = pathlib.Path('phylotype-set')
    mappings_dir = root / 'mappings'
    mappings_dir.mkdir(parents=True)
    with open('${placement_identity}') as in_h:
        identity = json.load(in_h)
    shutil.copyfile('${dedup_jplace}', root / 'cumulative.jplace')
    checksums = {}
    for source in sorted(pathlib.Path('.').glob('phylotype.*.parquet')):
        target = mappings_dir / source.name
        shutil.copyfile(source, target)
        checksums[str(target.relative_to(root))] = hashlib.sha256(target.read_bytes()).hexdigest()
    checksums['cumulative.jplace'] = hashlib.sha256((root / 'cumulative.jplace').read_bytes()).hexdigest()
    manifest = {
        'parent_set_id': None,
        'refpkg_id': identity['refpkg_id'],
        'tree_sha256': identity['tree_sha256'],
        'dataset_id': '${dataset_id}',
        'thresholds': ${thresholds},
        'distance': '${distance}',
        'lwr_overlap': ${lwr_overlap},
        'artifacts': checksums,
    }
    canonical = json.dumps(manifest, sort_keys=True, separators=(',', ':')).encode()
    manifest['set_id'] = hashlib.sha256(canonical).hexdigest()
    (root / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\\n')
    with tarfile.open('phylotype-set.tar.gz', 'w:gz') as out_h:
        out_h.add(root, arcname='phylotype-set')
    PY
    """
}
