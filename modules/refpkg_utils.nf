nextflow.enable.dsl=2

process ExtractRefpkg {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path refpkg_tgz

    output:
    path 'refpkg_tree.nwk', emit: tree
    path 'refpkg.aln.fasta', emit: ref_aln_fasta
    path 'refpkg.aln.sto', emit: ref_aln_sto
    path 'model.txt', emit: model
    path 'leaf_info.csv', emit: leaf_info
    path 'taxonomy.csv', emit: taxonomy
    path 'refpkg.cm', emit: cm
    path 'sv_registry.parquet', emit: sv_registry

    script:
    """
    python3 - <<'PY'
    import json
    import os
    import posixpath
    import re
    import tarfile

    from Bio import AlignIO

    archive = '${refpkg_tgz}'

    def fail(message):
        raise RuntimeError(f'Invalid refpkg: {message}')

    with tarfile.open(archive, 'r:*') as tar_h:
        members = {}
        for member in tar_h.getmembers():
            if not member.isfile():
                continue
            normalized = posixpath.normpath(member.name)
            if normalized.startswith('../') or normalized.startswith('/') or '/..' in normalized:
                fail(f'unsafe member {member.name!r}')
            basename = posixpath.basename(normalized)
            if basename in members:
                fail(f'duplicate member basename {basename!r}')
            members[basename] = member

        contents_member = members.get('CONTENTS.json')
        if contents_member is None:
            fail('missing CONTENTS.json')
        contents = json.load(tar_h.extractfile(contents_member))
        files = contents.get('files')
        if not isinstance(files, dict):
            fail('CONTENTS.json files is not an object')

        def extract(role):
            logical_name = files.get(role)
            if not logical_name:
                fail(f'missing logical component {role!r}')
            member = members.get(posixpath.basename(logical_name))
            if member is None:
                fail(f'component {role!r} is absent from archive')
            return tar_h.extractfile(member).read()

        registry_member = members.get('sv_registry.parquet')
        if registry_member is None:
            fail('missing sv_registry.parquet')

        for role in ('profile', 'tree', 'seq_info', 'taxonomy'):
            extract(role)

        with open('refpkg.cm', 'wb') as out_h:
            out_h.write(extract('profile'))
        with open('refpkg_tree.nwk', 'wb') as out_h:
            out_h.write(extract('tree'))
        with open('leaf_info.csv', 'wb') as out_h:
            out_h.write(extract('seq_info'))
        with open('taxonomy.csv', 'wb') as out_h:
            out_h.write(extract('taxonomy'))
        with open('sv_registry.parquet', 'wb') as out_h:
            out_h.write(tar_h.extractfile(registry_member).read())

        have_fasta = bool(files.get('aln_fasta'))
        have_sto = bool(files.get('aln_sto'))
        if not (have_fasta or have_sto):
            fail('missing both aln_fasta and aln_sto')
        if have_fasta:
            with open('refpkg.aln.fasta', 'wb') as out_h:
                out_h.write(extract('aln_fasta'))
        if have_sto:
            with open('refpkg.aln.sto', 'wb') as out_h:
                out_h.write(extract('aln_sto'))
        if not have_sto:
            with open('refpkg.aln.fasta', 'rt') as in_h, open('refpkg.aln.sto', 'wt') as out_h:
                AlignIO.write(AlignIO.read(in_h, 'fasta'), out_h, 'stockholm')
        if not have_fasta:
            with open('refpkg.aln.sto', 'rt') as in_h, open('refpkg.aln.fasta', 'wt') as out_h:
                AlignIO.write(AlignIO.read(in_h, 'stockholm'), out_h, 'fasta')

        if files.get('raxml_ng_model'):
            model = extract('raxml_ng_model').decode('utf-8')
        else:
            phylo_model = json.loads(extract('phylo_model').decode('utf-8'))
            rates = phylo_model.get('subs_rates', phylo_model)
            stats = extract('tree_stats').decode('utf-8')
            base_freqs = re.search(
                r'Base frequencies: (?P<A>0\\.\\d+) (?P<C>0\\.\\d+) (?P<G>0\\.\\d+) (?P<T>0\\.\\d+)',
                stats,
            )
            if base_freqs is None:
                fail('could not parse base frequencies from tree_stats')
            model = 'GTR{{{}}}+FU{{{}}}'.format(
                '/'.join(str(rates[key]) for key in ('ac', 'ag', 'at', 'cg', 'ct', 'gt')),
                '/'.join(base_freqs.group(key) for key in ('A', 'C', 'G', 'T')),
            )
        with open('model.txt', 'wt') as out_h:
            out_h.write(model)
    PY
    """
}
