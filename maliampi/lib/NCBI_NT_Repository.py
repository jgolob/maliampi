""" A wrapper class around either a database or
    filesystem-based cache of NCBI_NT.
"""
import logging
# Oh python 2.x vs 3.x
try:
    from urlparse import urlsplit, urljoin
except ImportError:
    from urllib.parse import urlsplit, urljoin

import pymongo
import gzip
import uuid

NT_REPO_NS = uuid.UUID('a3c2d868-1544-40b1-b341-278dcb756326')


class NCBI_NT_Repository(object):
    __engine__ = None  # Engine for the store. Right now only mongodb is supported
    __store_url__ = None

    COMPLEMENT = {
        'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
        'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
    }

    def reverse_complement(self, seq):
        bases = list(seq)
        bases = reversed([self.COMPLEMENT.get(base, base) for base in bases])
        bases = ''.join(bases)
        return bases

    def __init_mongodb__(self, ncbi_repository_db='ncbi_repo', collection='records'):
        import pymongo
        # Open the client
        mongodb_client = pymongo.MongoClient(self.__store_url__)
        self.__mongodb_client__ = mongodb_client

        # And our database, using the parameter above
        # (this will silently create the database if it doesn't exist)
        db = mongodb_client[ncbi_repository_db]
        self.__mongodb_db__ = db

        # And our collection
        collection = db[collection]
        self.__mongodb_collection__ = collection

        collection.create_index(
            'version',
            unique=True
        )

        # Peptides
        peptides = db['peptides']
        self.__mongodb_peptides__ = peptides

        peptides.create_index(
            [('seq', pymongo.HASHED)],
        )

        # Peptides
        rRNA16S = db['rRNA16S']
        self.__mongodb_rRNA16S__ = rRNA16S
        rRNA16S.create_index(
            [('seq', pymongo.HASHED)],
        )

        # Each record within the collection has a basic structure:
        """
            {'_id': ObjectId('5a4ff8c22d49d3afa1103e44'),
            # Every entry will have at least these three
            'accession': 'NC_000962',
            'version': 'NC_000962.3',
            'version_num': '3',

            # Some flags to help with quick filtering
            'is_genome': True,  # Is this accession a genome
            'genome_complete_checkm': 90  # Genome checkm completeness score
            'is_type': False,
            'taxonomy_verified': True, # Is the taxonomy verified by deenurp, etc

            # Most SHOULD have these eventually
            'tax_id': '83332',
            'feature_table': [(start, end, feature, qual, qual_val), ]
            'name': 'NC_000962',
            'organism': 'Mycobacterium tuberculosis H37Rv',
            'length': 1537,
            'modified_date': '14-DEC-2017',
            'mol_type': 'genomic DNA',
            'source': 'Mycobacterium tuberculosis H37Rv',
            'strain': 'H37Rv',
            'keywords': 'RefSeq;complete genome',
            'description': 'Mycobacterium tuberculosis H37Rv, complete genome',
            'isolate': '',
            'isolation_source': '',
            'pubmed_refs': ['20980199', '12368430', '9634230'],
            'refseq__accession': 'AL123456',
            'refseq__gi': '57116681',
            'refseq__seq_start': None,
            'refseq__seq_stop': None,

            # SOME will have (gzipped)
            'seq': 'TTTTGTTTG',
            }
        """
        # The very bare minimum is an accession, version, and version num.
        # tax_id is also frequently used.
        # We can create indicies for each of these.
        # version should be unique. Only one entry per version
        collection.create_index('version', unique=True)
        # We will be frequently selecting all of a given tax ID
        collection.create_index('tax_id', background=True)

    def __mongodb_get_existing_versions__(self, contraints):
        collection = self.__mongodb_collection__
        try:
            return set(collection.distinct('version'))
        except pymongo.errors.OperationFailure:
            return {r['version'] for r in collection.find(contraints, {'version': 1, '_id': 0})}

    def __mongodb_add_entry__(self, entry):
        collection = self.__mongodb_collection__
        try:
            collection.update_one(
                filter={'version': entry['version']},
                update={'$set': entry},
                upsert=True
            )
        except pymongo.errors.DocumentTooLarge:
            logging.warning("Failed to save document for {} due to size. Reducing".format(
                entry['version']))
            entry['feature_table'] = [
                f for f in entry.get('feature_table', [])
                if f.get('feature') == 'rRNA'
            ]
            try:
                collection.update_one(
                    filter={'version': entry['version']},
                    update={'$set': entry},
                    upsert=True
                )
            except pymongo.errors.DocumentTooLarge:
                logging.error("Failed to save document for {}  after reducing".format(
                    entry['version']))

    def __mongodb_versions_needing_feature_table__(self, constraints={}):
        collection = self.__mongodb_collection__
        constraints.update(
            {'feature_table': None}
        )
        return {
            r['version'] for r in
            collection.find(constraints, {'version': 1, "_id": 0})
        }

    def __mongodb_versions_needing_data__(self, constraints={}):
        collection = self.__mongodb_collection__
        constraints.update(
            {'modified_date': None}
        )
        return {
            r['version'] for r in
            collection.find(constraints, {'version': 1, "_id": 0})
        }

    def __mongodb_versions_needing_peptides__(self, constraints={}):
        collection = self.__mongodb_collection__
        constraints.update({
            'peptides': None,
            'is_genome': True
            })
        return {
            r['version'] for r in
            collection.find(constraints, {'version': 1, "_id": 0})
        }

    def __mongodb_add_feature__(self, feature_dict):
        collection = self.__mongodb_collection__
        collection.update_one(
            {'version': feature_dict['id']},  # constraint
            {
                '$addToSet': {'feature_table': (
                    feature_dict.get('feature'),
                    feature_dict.get('qual'),
                    feature_dict.get('qual_value'),
                    feature_dict.get('seq_start'),
                    feature_dict.get('seq_stop'),
                    feature_dict.get('strand')
                    )}
            }
        )

    def __mongodb_update_feature_table__(self, version, feature_table):
        collection = self.__mongodb_collection__
        try:
            collection.update_one(
                {'version': version},  # constraint
                {
                    '$set': {'feature_table': feature_table}
                }
            )
        except pymongo.errors.DocumentTooLarge:
            logging.error(
                "Feature table makes the document too large. Limiting to rRNA, CDS, and gene"
                )
            feature_table = [f for f in feature_table if
                             (f.get('feature') == 'rRNA' or
                              f.get('feature') == 'CDS' or
                              f.get('feature') == 'gene')]
            try:
                collection.update_one(
                    {'version': version},  # constraint
                    {
                        '$set': {'feature_table': feature_table}
                    }
                )
            except:
                logging.error("Failed to update feature table on {}".format(version))

    def __mongodb_find_feature__(
            self,
            version,
            feature,
            qual,
            qual_val,
            extra_qualifiers={},
            extra_projection={}):
        collection = self.__mongodb_collection__
        qualifiers = {}
        ft_qualifiers = {}
        if version:
            qualifiers['version'] = version
        if feature:
            ft_qualifiers['feature'] = feature
        if qual_val:
            ft_qualifiers['qual_val'] = qual_val
        if qual:
            ft_qualifiers['qual'] = qual
        qualifiers['feature_table'] = {'$elemMatch': ft_qualifiers}
        projection = {'version': 1, 'feature_table': 1, 'seq': 1, '_id': 0}
        projection.update(extra_projection)

        for rec_match in collection.find(
                filter=qualifiers,
                projection=projection):
            rec_seq = gzip.decompress(rec_match['seq'])
            if feature:
                f_match = [f for f in rec_match['feature_table'] if f['feature'] == feature]
            else:
                f_match = rec_match['feature_table']
            if qual:
                q_match = [f for f in f_match if f['qual'] == qual]
            else:
                q_match = f_match
            if qual_val:
                qv_match = [f for f in q_match if f['qual_val'] == qual_val]
            else:
                qv_match = q_match

            for qv_m in qv_match:
                start_stop = sorted([qv_m['seq_start'], qv_m['seq_stop']])
                qv_seq = rec_seq[start_stop[0]:start_stop[1]].decode('utf-8')
                if qv_m['strand'] == -1:
                    qv_seq = self.reverse_complement(qv_seq)
                f_dict = {
                    'accession_version': rec_match['version'],
                    'start': start_stop[0],
                    'stop': start_stop[1],
                    'feature': qv_m['feature'],
                    'qual': qv_m['qual'],
                    'qual_val': qv_m['qual_val'],
                    'seq': qv_seq
                }
                f_dict.update({
                    k: rec_match.get(k)
                    for k in extra_projection.keys()
                })
                yield f_dict

    def __mongodb_get_full_sequence__(self, version):
        collection = self.__mongodb_collection__
        rec = collection.find_one({
            'version': version
        })
        if rec is None:
            logging.error("Could not find {}".format(version))
            return None
        # Implicit else
        seq = gzip.decompress(rec['seq']).decode('utf-8')
        if '>' not in seq:
            seq = ">{}\n".format(version)+seq
        return seq

    def __mongodb_add_peptides_to_version__(self, version, peptides):
        collection = self.__mongodb_collection__
        peptides_collection = self.__mongodb_peptides__

        collection.update_one(
            filter={'version': version},
            update={'$addToSet': {'peptides': {'$each': peptides}}},
        )
        for p in peptides:
            peptides_collection.update_one(
                filter={'seq': p},
                update={'$addToSet': {'genomes': version}},
                upsert=True
            )
            peptides_collection.update_one(
                filter={'seq': p},
                update={'$set': {'seq_id': uuid.uuid3(NT_REPO_NS, p)}},
                upsert=True
            )

    def __mongodb_add_rRNA16s_to_version__(self, version, rRNA16S):
        collection = self.__mongodb_collection__
        rRNA_collection = self.__mongodb_rRNA16S__
        collection.update_one(
            filter={'version': version},
            update={'$addToSet': {'rRNA_16S': {'$each': rRNA16S}}},
        )
        for r in rRNA16S:
            rRNA_collection.update_one(
                filter={'seq': r},
                update={'$addToSet': {'genomes': version}},
                upsert=True
            )
            rRNA_collection.update_one(
                filter={'seq': r},
                update={'$set': {'seq_id': uuid.uuid3(NT_REPO_NS, r)}},
                upsert=True
            )

    def get_full_sequence(self, version):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_get_full_sequence__(version)

    def find_feature(
            self,
            version=None,
            feature=None,
            qual=None,
            qual_val=None,
            extra_qualifiers={},
            extra_fields=[]):
        if not (version or feature or qual or qual_val):
            logging.warn("Nothing specified to find")
            return([])
        # Implicit else
        if self.__engine__ == 'mongodb':
            extra_projection = {f: 1 for f in extra_fields}
            return self.__mongodb_find_feature__(
                version=version,
                feature=feature,
                qual=qual,
                qual_val=qual_val,
                extra_qualifiers=extra_qualifiers,
                extra_projection=extra_projection)

    def update_feature_table(self, version, feature_table):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_update_feature_table__(version, feature_table)

    def add_feature(self, feature_dict):
        # Takes a dict with the minimum entries:
        # id: accession_version from which this feature came
        # seq_start: seq_start
        # seq_stop: seq_stop
        # strand: 1 or 2
        # feature: gene / cds etc
        # qual: product, gene, etc
        # qual_value: dnaA etc,
        #
        # and makes entries in that records feature table
        if self.__engine__ == 'mongodb':
            return self.__mongodb_add_feature__(feature_dict)

    def versions_needing_peptides(self, constraints={}):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_versions_needing_peptides__(constraints)

    def versions_needing_data(self, constraints={}):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_versions_needing_data__(constraints)

    def has_versions(self, versions):
        # Method that returns True if all versions are within repo. Else false
        return len(set(versions) - self.get_existing_versions()) == 0

    def versions_needing_feature_table(self, constraints={}):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_versions_needing_feature_table__(constraints)

    def add_entry(self, entry):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_add_entry__(entry)

    def add_peptides_to_version(self, version, peptides):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_add_peptides_to_version__(
                version,
                peptides
            )

    def add_rRNA16s_to_version(self, version, rRNA16S):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_add_rRNA16s_to_version__(
                version,
                rRNA16S
            )

    def get_existing_versions(self, contraints={}):
        if self.__engine__ == 'mongodb':
            return self.__mongodb_get_existing_versions__(contraints)

    def add_new_versions(self, input_versions, extra_values={}):
        input_versions_set = set(input_versions)
        existing_versions = self.get_existing_versions()
        new_versions = (input_versions_set - existing_versions)
        for ver in new_versions:
            entry = {
                'version': ver,
                'accession': ver.split('.')[0],
                'version_num': ver.split('.')[1],
            }
            entry.update(extra_values)
            self.add_entry(entry)
        # For convenience, return the new versions
        return new_versions

    def __init__(self, store_url, ncbi_repository_db='ncbi_repo', collection='records'):

        self.__store_url__ = store_url
        # Figure out where our store is, using urlsplit
        store_split = urlsplit(store_url)
        if store_split.scheme == 'mongodb':
            self.__engine__ = 'mongodb'
            self.__init_mongodb__(ncbi_repository_db, collection)
        else:
            raise Exception("Engine {} is not supported currently".format(
                store_split.scheme
            ))
