"""
Class to Extract Genbank records and return as a dict
"""

import re
from datetime import date
from Bio import SeqIO
import logging
import gzip


class ExtractGenbank(object):
    # columns of output files
    ANNOTATION_COLS = ['seqname', 'version', 'accession', 'name',
                       'description', 'tax_id', 'modified_date', 'download_date',
                       'version_num', 'source', 'keywords', 'organism', 'length',
                       'ambig_count', 'strain', 'mol_type', 'isolate',
                       'isolation_source', 'seq_start', 'seq_stop']

    PUBMED_COLS = ['pubmed_id', 'version', 'accession']

    REFERENCE_COLS = ['pubmed_id', 'title', 'authors',
                      'journal', 'consrtm', 'comment']

    REFSEQ_INFO_COLS = ['seqname', 'accession', 'gi', 'seq_start', 'seq_stop']

    REFSEQ_SOURCE = re.compile(r'(?P<accession>[A-Z]{1,4}\d{5,8})'
                               r':?(?P<seq_start>\d+)?-?(?P<seq_stop>\d+)?')

    GI_SOURCE = re.compile(r'gi:(?P<gi>\d+)')

    ACGT = frozenset('ACGT')

    # Key is version. Value is record in dict format
    records = {}

    def __init__(self,
                 record_file,
                 gzip_seq=True,
                 focused_feature_table=True,
                 download_date=str(date.today())):

        download_date = download_date
        try:
            raw_records = SeqIO.parse(record_file, 'genbank')

            for rec_i, raw_rec in enumerate(raw_records):
                try:
                    record = self.parse_record(
                        raw_rec,
                        focused_feature_table=focused_feature_table)

                    record.update({'download_date': download_date})
                    if gzip_seq:
                        record.update({'seq': gzip.compress(
                            str(raw_rec.seq).encode('utf8')
                        )})
                    else:
                        record.update({'seq': str(raw_rec.seq)})
                    refs = self.parse_references(raw_rec)
                    record.update({'pubmed_refs': [ref['pubmed_id'] for ref in refs]})

                    # refseqs
                    if '_' in record['accession']:
                        refseq = self.parse_refseq_source(raw_rec)
                        record.update({"refseq__{}".format(rsk): refseq[rsk] for rsk in refseq})

                    self.records[record['version']] = record

                except Exception as e:
                    logging.error("{} when working on {}th record".format(
                        e,
                        rec_i,
                    ))
        except Exception as e:
            logging.error("Failed to load {} due to error {}".format(record_file, e))
            raw_records = []

    def get_records(self):
        return self.records

    def parse_coordinates(self, g):
        accessions = g.annotations['accessions']
        if 'REGION:' in accessions:
            coordinates = accessions[accessions.index('REGION:') + 1]
            matches = re.findall(r'\d+', coordinates)
            if len(matches) == 1:
                # some records are strange...
                seq_start, seq_stop = matches[0], len(g.seq)
            else:
                seq_start, seq_stop = matches
        else:
            seq_start, seq_stop = 1, len(g.seq)
        return seq_start, seq_stop

    def parse_record(self, g, focused_feature_table):
        version, accession, version_num = self.parse_version(g)

        # features
        if focused_feature_table:
            # Make a feature table only covering CDS, genes, and rRNA
            # output is a LIST where each entry is a dict
            # 'seq_start', 'seq_stop', 'strand', 'feature', 'qual', 'qual_value' at a min
            # Additional values are stored
            feature_table = []
            for feature in (f for f in g.features
                            if (f.type == 'CDS' or f.type == 'rRNA')):
                feature_dict = {
                    'seq_start': int(feature.location.start),
                    'seq_stop': int(feature.location.end),
                    'strand': int(feature.location.strand),
                    'feature': feature.type
                }
                qual_dict = dict(feature.qualifiers)
                try:
                    del qual_dict['translation']
                except KeyError:
                    pass
                for qual_key, qual_value in qual_dict.items():
                    for val in qual_value:
                        qdr = feature_dict.copy()
                        qdr.update({'qual': qual_key,
                                    'qual_val': val})
                        feature_table.append(qdr)

        source = next(i for i in g.features if i.type == 'source')
        quals = source.qualifiers

        # get tax_id
        if 'db_xref' in quals:
            for i in quals.get('db_xref', []):
                if i.startswith('taxon:'):
                    tax_id = i[6:]
                    break
        # occasionally, tax_ids are missing
        else:
            tax_id = ''

        seq_start, seq_stop = self.parse_coordinates(g)
        if all([seq_start, seq_stop]):
            seq_id = '{}_{}_{}'.format(accession, seq_start, seq_stop)
        else:
            seq_id = g.id

        info = dict(accession=accession,
                    ambig_count=sum(1 for b in g.seq if b not in self.ACGT),
                    modified_date=g.annotations['date'],
                    description=g.description,
                    keywords=';'.join(g.annotations.get('keywords', [])),
                    length=len(g),
                    name=g.name,
                    organism=g.annotations['organism'],
                    seq_start=seq_start,
                    seq_stop=seq_stop,
                    seqname=seq_id,
                    source=g.annotations['source'],
                    mol_type=';'.join(quals.get('mol_type', '')),
                    strain=';'.join(quals.get('strain', '')),
                    isolate=';'.join(quals.get('isolate', '')),
                    isolation_source=';'.join(quals.get('isolation_source', '')),
                    tax_id=tax_id,
                    version=version,
                    version_num=version_num)
        if focused_feature_table:
            info['feature_table'] = feature_table

        return info

    def parse_references(self, g):
        """
        Parse reference annotations that have a pubmed_id
        """
        references = []
        if 'references' in g.annotations:
            refs = [r for r in g.annotations['references'] if r.pubmed_id]
            for r in refs:
                references.append(
                    dict(title=r.title,
                         authors=r.authors,
                         comment=r.comment,
                         consrtm=r.consrtm,
                         journal=r.journal,
                         pubmed_id=r.pubmed_id))
        return references

    def parse_refseq_source(self, g):
        acc = re.search(self.REFSEQ_SOURCE, g.annotations['comment'])
        gi = re.search(self.GI_SOURCE, g.annotations['comment'])

        if not acc and not gi:
            raise ValueError('Cannot parse record')

        result = {}
        if acc:
            result.update(acc.groupdict())
        if gi:
            result.update(gi.groupdict())

        return result

    def parse_version(self, g):
        """
        Return the accession and version of a Bio.SeqRecord.SeqRecord
        """
        annotations = g.annotations
        accession = annotations.get('accessions', [''])[0]
        if accession:
            if 'sequence_version' in annotations:
                version_num = str(annotations.get('sequence_version'))
            elif g.id.startswith(accession + '.'):
                version_num = g.id.split('.', 1)[1]
            else:
                version_num = '1'
            version = accession + '.' + version_num
        else:
            version = ''
        return version, accession, version_num
