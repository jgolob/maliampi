nextflow.enable.dsl=2

/* Canonicalize the one abundance representation used by downstream stages. */

process NormalizeSvLong {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'
    publishDir "${params.output}/sv", mode: 'copy'

    input:
    path sv_fasta
    path sv_long
    val dataset_id

    output:
    path 'sv.h5ad', emit: sv_h5ad

    script:
    """
    python3 - <<'PY'
    import csv
    import re

    dataset_id = '${dataset_id}'
    if dataset_id and not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', dataset_id):
        raise SystemExit('dataset_id must be a non-empty safe identifier')

    def namespaced(value):
        return f'{dataset_id}__{value}' if dataset_id else value

    records = []
    seen = set()
    header = None
    seq = []
    with open('${sv_fasta}') as in_h:
        for line in in_h:
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header, seq = line[1:].strip().split()[0], []
                if header in seen:
                    raise SystemExit(f'duplicate FASTA ID: {header}')
                seen.add(header)
            else:
                seq.append(line.strip())
        if header is not None:
            records.append((header, ''.join(seq)))
    if not records:
        raise SystemExit('SV FASTA is empty')

    rows = []
    with open('${sv_long}', newline='') as in_h:
        reader = csv.DictReader(in_h)
        if not {'specimen', 'sv', 'count'} <= set(reader.fieldnames or []):
            raise SystemExit('sv_long requires specimen, sv, count columns')
        for row in reader:
            if not row['specimen'] or not row['sv']:
                raise SystemExit('sv_long has an empty specimen or SV identifier')
            try:
                count = int(row['count'])
            except ValueError:
                raise SystemExit(f"invalid count for {row['sv']}")
            if count < 0:
                raise SystemExit(f"negative count for {row['sv']}")
            rows.append((row['specimen'], row['sv'], count))
    long_ids = {sv for _, sv, count in rows if count > 0}
    if long_ids != seen:
        missing = sorted(seen - long_ids)
        extra = sorted(long_ids - seen)
        raise SystemExit(f'FASTA/SV-long mismatch; missing={missing}, extra={extra}')

    with open('sv.fasta', 'w') as out_h:
        for sv, sequence in records:
            out_h.write(f'>{namespaced(sv)}\\n{sequence}\\n')
    with open('sv.long.csv', 'w', newline='') as out_h:
        writer = csv.writer(out_h)
        writer.writerow(('specimen', 'sv', 'count'))
        for specimen, sv, count in sorted(rows):
            writer.writerow((specimen, namespaced(sv), count))

    # pplacer's historic adapter requires one representative specimen-SV ID per
    # specimen/SV row. Stable ordering makes this reproducible.
    representative = {}
    for specimen, sv, count in sorted(rows, key=lambda r: (-r[2], r[0], r[1])):
        representative.setdefault(sv, specimen)
    with open('sv.map.csv', 'w', newline='') as map_h, open('sv.weights.csv', 'w', newline='') as weights_h:
        map_writer, weight_writer = csv.writer(map_h), csv.writer(weights_h)
        for specimen, sv, count in sorted(rows):
            sv_sp = sv if representative[sv] == specimen else f'{sv}__{specimen}'
            sv_sp = namespaced(sv_sp)
            map_writer.writerow((sv_sp, specimen))
            weight_writer.writerow((namespaced(sv), sv_sp, count))
    PY
    maliampi-sv-import --fasta ${sv_fasta} --long ${sv_long} \
      --project-id ${dataset_id ?: 'legacy'} --dataset-id ${dataset_id ?: 'legacy'} \
      --output-h5ad sv.h5ad
    """
}

process PplacerReduplicate {
    container "${params.container__pplacer}"
    label 'io_limited'
    publishDir "${params.output}/placement", mode: 'copy'

    input:
    path dedup_jplace
    path sv_weights

    output:
    path 'redup.jplace.gz', emit: redup_jplace

    script:
    """
    guppy redup -m -o /dev/stdout -d ${sv_weights} ${dedup_jplace} | gzip > redup.jplace.gz
    """
}

process SharetableToSvLong {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path sharetable

    output:
    path 'sharetable.sv.long.csv', emit: long

    script:
    """
    python3 - <<'PY'
    import csv

    with open('${sharetable}', newline='') as in_h:
        reader = csv.DictReader(in_h, delimiter='\\t')
        fields = reader.fieldnames or []
        reserved = {'label', 'group', 'numotus', 'numsvs'}
        variants = [field for field in fields if field.lower() not in reserved]
        if not variants or 'group' not in {field.lower() for field in fields}:
            raise SystemExit('sharetable must contain group and one or more SV columns')
        group_field = next(field for field in fields if field.lower() == 'group')
        with open('sharetable.sv.long.csv', 'w', newline='') as out_h:
            writer = csv.writer(out_h)
            writer.writerow(('specimen', 'sv', 'count'))
            for row in reader:
                specimen = row[group_field]
                for sv in variants:
                    try:
                        count = int(row[sv])
                    except ValueError:
                        raise SystemExit(f'invalid count for {specimen}/{sv}')
                    if count < 0:
                        raise SystemExit(f'negative count for {specimen}/{sv}')
                    if count:
                        writer.writerow((specimen, sv, count))
    PY
    """
}

process MapWeightsToSvLong {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path sv_map
    path sv_weights

    output:
    path 'map_weights.sv.long.csv', emit: long

    script:
    """
    python3 - <<'PY'
    import csv

    mapping = {}
    with open('${sv_map}', newline='') as in_h:
        for row in csv.reader(in_h):
            if len(row) != 2 or not all(row):
                raise SystemExit('sv_map must be headerless (specimen_sv_id,specimen)')
            if row[0] in mapping and mapping[row[0]] != row[1]:
                raise SystemExit(f'conflicting specimen assignment for {row[0]}')
            mapping[row[0]] = row[1]
    counts = {}
    with open('${sv_weights}', newline='') as in_h:
        for row in csv.reader(in_h):
            if len(row) != 3 or not all(row):
                raise SystemExit('sv_weights must be headerless (sv,specimen_sv_id,count)')
            sv, specimen_sv, count = row
            if specimen_sv not in mapping:
                raise SystemExit(f'weight references unmapped specimen SV {specimen_sv}')
            try:
                count = int(count)
            except ValueError:
                raise SystemExit(f'invalid count for {sv}')
            if count < 0:
                raise SystemExit(f'negative count for {sv}')
            key = (mapping[specimen_sv], sv)
            counts[key] = counts.get(key, 0) + count
    with open('map_weights.sv.long.csv', 'w', newline='') as out_h:
        writer = csv.writer(out_h)
        writer.writerow(('specimen', 'sv', 'count'))
        writer.writerows((specimen, sv, count) for (specimen, sv), count in sorted(counts.items()) if count)
    PY
    """
}

process ValidateSvMapWeights {
    container "${params.container__maliampi_tools}"
    label 'maliampi_tools'
    label 'io_limited'

    input:
    path sv_fasta
    path sv_long
    path sv_map
    path sv_weights

    output:
    path 'checked.fasta', emit: fasta
    path 'checked.long.csv', emit: long

    script:
    """
    python3 - <<'PY'
    import csv
    import shutil

    expected = {}
    with open('${sv_long}', newline='') as in_h:
        reader = csv.DictReader(in_h)
        if not {'specimen', 'sv', 'count'} <= set(reader.fieldnames or []):
            raise SystemExit('sv_long requires specimen, sv, count columns')
        for row in reader:
            expected[(row['specimen'], row['sv'])] = expected.get((row['specimen'], row['sv']), 0) + int(row['count'])
    mapping = {row[0]: row[1] for row in csv.reader(open('${sv_map}', newline='')) if len(row) == 2}
    observed = {}
    for row in csv.reader(open('${sv_weights}', newline='')):
        if len(row) != 3 or row[1] not in mapping:
            raise SystemExit('weights/map inputs are not a valid canonical pair')
        key = (mapping[row[1]], row[0])
        observed[key] = observed.get(key, 0) + int(row[2])
    if expected != observed:
        raise SystemExit('sv_long does not agree with sv_map/sv_weights')
    shutil.copyfile('${sv_fasta}', 'checked.fasta')
    shutil.copyfile('${sv_long}', 'checked.long.csv')
    PY
    """
}
