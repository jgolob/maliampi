nextflow.enable.dsl=2

process ValidateRefpkg {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path refpkg_tgz

    output:
    path 'validated_refpkg.tar.gz', emit: refpkg
    path 'refpkg.identity.json', emit: identity
    path 'refpkg.validation.json', emit: report

    script:
    """
    python3 - <<'PY'
    import csv
    import hashlib
    import json
    import math
    import posixpath
    import re
    import shutil
    import tarfile

    archive = '${refpkg_tgz}'

    def fail(message):
        raise RuntimeError(f'Invalid refpkg: {message}')

    def sha256(data):
        return hashlib.sha256(data).hexdigest()

    def parse_fasta(data):
        records = {}
        current = None
        for line in data.decode('utf-8').splitlines():
            if line.startswith('>'):
                current = line[1:].split()[0]
                if not current or current in records:
                    fail('invalid or duplicate FASTA identifier')
                records[current] = []
            elif line.strip():
                if current is None:
                    fail('FASTA sequence before identifier')
                records[current].append(line.strip())
        if not records:
            fail('empty FASTA alignment')
        return {key: ''.join(value) for key, value in records.items()}

    def parse_stockholm(data):
        records = {}
        for line in data.decode('utf-8').splitlines():
            line = line.strip()
            if not line or line.startswith('#') or line == '//':
                continue
            fields = line.split()
            if len(fields) != 2:
                fail('invalid Stockholm alignment row')
            records.setdefault(fields[0], []).append(fields[1])
        if not records:
            fail('empty Stockholm alignment')
        return {key: ''.join(value) for key, value in records.items()}

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
        expected_md5 = contents.get('md5', {})

        def component(role):
            logical_name = files.get(role)
            if not logical_name:
                fail(f'missing logical component {role!r}')
            member = members.get(posixpath.basename(logical_name))
            if member is None:
                fail(f'component {role!r} is absent from archive')
            data = tar_h.extractfile(member).read()
            expected = expected_md5.get(role)
            if expected and hashlib.md5(data).hexdigest() != expected:
                fail(f'MD5 mismatch for {role!r}')
            return data

        required = ('profile', 'tree', 'seq_info', 'taxonomy')
        components = {role: component(role) for role in required}
        registry_member = members.get('sv_registry.parquet')
        if registry_member is None:
            fail('missing sv_registry.parquet')
        components['sv_registry'] = tar_h.extractfile(registry_member).read()
        if files.get('aln_fasta'):
            components['aln_fasta'] = component('aln_fasta')
        if files.get('aln_sto'):
            components['aln_sto'] = component('aln_sto')
        if not {'aln_fasta', 'aln_sto'} & components.keys():
            fail('missing both aln_fasta and aln_sto')
        if files.get('raxml_ng_model'):
            components['raxml_ng_model'] = component('raxml_ng_model')
        else:
            components['phylo_model'] = component('phylo_model')
            components['tree_stats'] = component('tree_stats')

        fasta = parse_fasta(components['aln_fasta']) if 'aln_fasta' in components else None
        sto = parse_stockholm(components['aln_sto']) if 'aln_sto' in components else None
        alignment = fasta or sto
        if fasta and sto:
            if set(fasta) != set(sto):
                fail('FASTA and Stockholm identifiers differ')
            if {key: len(value) for key, value in fasta.items()} != {key: len(value) for key, value in sto.items()}:
                fail('FASTA and Stockholm alignment lengths differ')

        tree_text = components['tree'].decode('utf-8').strip()
        leaves = set(re.findall(r'(?<=[(,])([^():,;{}\\s]+)(?=:)', tree_text))
        if not leaves:
            fail('could not parse tree leaves')
        if leaves != set(alignment):
            fail('tree leaves and alignment identifiers differ')

        seq_info = list(csv.DictReader(components['seq_info'].decode('utf-8').splitlines()))
        seq_names = {row.get('seqname') for row in seq_info}
        if leaves - seq_names:
            fail('sequence information does not cover tree leaves')
        tax_ids = {row.get('tax_id') for row in seq_info if row.get('seqname') in leaves}
        taxonomy = list(csv.DictReader(components['taxonomy'].decode('utf-8').splitlines()))
        taxonomy_ids = {row.get('tax_id') for row in taxonomy}
        if tax_ids - taxonomy_ids:
            fail('taxonomy does not cover sequence-information tax IDs')

        if not components['profile'].strip():
            fail('empty covariance model')
        model_role = 'raxml_ng_model' if 'raxml_ng_model' in components else 'phylo_model'
        if not components[model_role].strip():
            fail('empty placement model')

        digest_rows = [f'{role}:{sha256(data)}' for role, data in sorted(components.items())]
        refpkg_id = sha256(('\\n'.join(digest_rows) + '\\n').encode('utf-8'))
        identity = {
            'format_version': 1,
            'refpkg_id': refpkg_id,
            'components': {role: sha256(data) for role, data in sorted(components.items())},
        }
        report = {
            'valid': True,
            'refpkg_id': refpkg_id,
            'tree_leaf_count': len(leaves),
            'alignment_length': len(next(iter(alignment.values()))),
            'metadata': contents.get('metadata', {}),
        }

    shutil.copyfile(archive, 'validated_refpkg.tar.gz')
    with open('refpkg.identity.json', 'wt') as out_h:
        json.dump(identity, out_h, indent=2, sort_keys=True)
        out_h.write('\\n')
    with open('refpkg.validation.json', 'wt') as out_h:
        json.dump(report, out_h, indent=2, sort_keys=True)
        out_h.write('\\n')
    PY
    """
}
