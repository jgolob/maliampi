import sciluigi as sl
import luigi
from .NCBI_NT_Repository import NCBI_NT_Repository
from tempfile import NamedTemporaryFile
import sqlite3
import logging
import tarfile
import io
import os
import json
import csv
import shlex
import argparse

log = logging.getLogger('sciluigi-interface')


# BASE CLASS #


class NCBI_Repo_TargetInfo(sl.TargetInfo):
    task = None
    target = None
    path = None

    def __init__(self, task, path):
        self.task = task
        self.path = path
        self.target = NCBI_Repo_Target
        self.target.path = path
        self.target.repo = NCBI_NT_Repository(path)

    def open(self):
        raise NotImplementedError("Open is not implemented for NCBI_Repo_TargetInfo")


class NCBI_Repo_Target(luigi.Target):
    # A class to wrap around NCBI_NT_Repository
    path = None
    repo = None

    @classmethod
    def exists(self):
        # Simple. So long as init succeded return True. Not very useful.
        return self.repo

###
# Check to see we have every entry in a given list.


class NCBI_Repo_Entries_TargetInfo(sl.TargetInfo):
    def __init__(self, task, path, all_entries):
        self.task = task
        self.path = path
        self.target = NCBI_Repo_Entries_Target
        self.target.repo = NCBI_NT_Repository(path)
        self.target.all_entries = all_entries

    def open(self):
        raise NotImplementedError("Open is not implemented for NCBI_Repo_TargetInfo")


class NCBI_Repo_Entries_Target(luigi.Target):
    # A class to wrap around NCBI_NT_Repository and see if all the
    # accession versions are already included
    repo = None
    all_entries = None

    @classmethod
    def exists(self):
        if self.all_entries == -1:
            return False
        # Make a set of all_versions. Get the existing versions from the repo.
        # Find out the difference. If the length is zero, we have everything.
        return self.repo.has_versions(self.all_entries)

###
# Check if every entry has been downloaded


class NCBI_Repo_Filled_TargetInfo(sl.TargetInfo):
    def __init__(self, task, path):
        self.task = task
        self.path = path
        self.target = NCBI_Repo_Filled_Target
        self.target.repo = NCBI_NT_Repository(path)

    def open(self):
        raise NotImplementedError("Open is not implemented for NCBI_Repo_TargetInfo")


class NCBI_Repo_Filled_Target(luigi.Target):
    # A class to wrap around NCBI_NT_Repository and see if all the
    # every entry has a FT
    repo = None

    @classmethod
    def exists(self):
        # Make a set of all_versions. Get the existing versions from the repo.
        # Find out the difference. If the length is zero, we have everything.
        return (len(self.repo.versions_needing_data()) == 0)


###
# Check if every genome entry has had peptides and rRNA 16S extracted

class NCBI_Repo_Peptides_TargetInfo(sl.TargetInfo):
    def __init__(self, task, path):
        self.task = task
        self.path = path
        self.target = NCBI_Repo_Peptides_Target
        self.target.repo = NCBI_NT_Repository(path)

    def open(self):
        raise NotImplementedError("Open is not implemented for NCBI_Repo_TargetInfo")


class NCBI_Repo_Peptides_Target(luigi.Target):
    # A class to wrap around NCBI_NT_Repository and see if all the
    # every entry has a FT
    repo = None

    @classmethod
    def exists(self):
        # Make a set of all_versions. Get the existing versions from the repo.
        # Find out the difference. If the length is zero, we have everything.
        return (len(self.repo.versions_needing_peptides()) == 0)


# Refpkg.tgz
class RefpkgTGZ_ContainerTargetInfo(sl.ContainerTargetInfo):
    # A light wrapper to support a tar.gz refpkg
    def __init__(self, *args, **kwargs):
        super(RefpkgTGZ_ContainerTargetInfo, self).__init__(
            *args,
            **kwargs
        )
        self.format = luigi.format.Gzip
        self.target.format = luigi.format.Gzip

    def get_aln_loc_fasta(self):
        with self.target.open('rb') as refpkg_tgz_h:
            tar_h = tarfile.open(
                mode='r:*',
                fileobj=io.BytesIO(refpkg_tgz_h.read())
            )
            # Get the contents keyed by the filenames, excluding dirs
            tar_contents_dict = {os.path.basename(f.name): f for f in tar_h.getmembers()}
            contents = json.load(tar_h.extractfile(tar_contents_dict['CONTENTS.json']))
            return contents['files'].get('aln_fasta')

    def get_aln_loc_sto(self):
        with self.target.open('rb') as refpkg_tgz_h:
            tar_h = tarfile.open(
                mode='r:*',
                fileobj=io.BytesIO(refpkg_tgz_h.read())
            )
            # Get the contents keyed by the filenames, excluding dirs
            tar_contents_dict = {os.path.basename(f.name): f for f in tar_h.getmembers()}
            contents = json.load(tar_h.extractfile(tar_contents_dict['CONTENTS.json']))
            return contents['files'].get('aln_sto')

    def get_refpkg_rel_path(self):
        with self.target.open('rb') as refpkg_tgz_h:
            refpkg_tar_h = tarfile.open(
                mode='r:*',
                fileobj=io.BytesIO(refpkg_tgz_h.read())
            )
            return os.path.dirname(refpkg_tar_h.getnames()[0])

    def get_tax_ids(self):
        with self.target.open('rb') as refpkg_tgz_h:
            tar_h = tarfile.open(
                mode='r:*',
                fileobj=io.BytesIO(refpkg_tgz_h.read())
            )
            # Get the contents keyed by the filenames, excluding dirs
            tar_contents_dict = {os.path.basename(f.name): f for f in tar_h.getmembers()}
            contents = json.load(tar_h.extractfile(tar_contents_dict['CONTENTS.json']))
            taxonomy_h = tar_h.extractfile(
                os.path.join(
                    self.get_refpkg_rel_path(),
                    contents['files']['taxonomy']
                )
            )
            taxonomy_reader = csv.DictReader(
                io.TextIOWrapper(taxonomy_h)
                )
            return {r.get('tax_id') for r in taxonomy_reader}

    def get_tax_ranks(self):
        with self.target.open('rb') as refpkg_tgz_h:
            tar_h = tarfile.open(
                mode='r:*',
                fileobj=io.BytesIO(refpkg_tgz_h.read())
            )
            # Get the contents keyed by the filenames, excluding dirs
            tar_contents_dict = {os.path.basename(f.name): f for f in tar_h.getmembers()}
            contents = json.load(tar_h.extractfile(tar_contents_dict['CONTENTS.json']))
            taxonomy_h = tar_h.extractfile(
                os.path.join(
                    self.get_refpkg_rel_path(),
                    contents['files']['taxonomy']
                )
            )
            taxonomy_reader = csv.DictReader(
                io.TextIOWrapper(taxonomy_h)
                )
            return {r.get('rank') for r in taxonomy_reader}

#  Placement DB


class PlacementDB_Prepped_ContainerTargetInfo(sl.ContainerTargetInfo):
    def __init__(
            self,
            task,
            path,
            is_tmp=False,
            client=None,
            removed_malformed=True,
            expected_ranks=None,
            expected_tax_ids=None,
            ):
                super(PlacementDB_Prepped_ContainerTargetInfo, self).__init__(
                    task=task,
                    path=path,
                    format=luigi.format.Nop,
                    is_tmp=is_tmp,
                    client=client
                    )
                self.removed_malformed = removed_malformed
                self.target.basic_exists = self.target.exists
                self.target.exists = self.full_exists
                self.expected_ranks = expected_ranks
                self.expected_tax_ids = expected_tax_ids

    def full_exists(self):
        if not self.target.basic_exists():
            return False
        # Implicit else, the file exists. Check our db for completeness
        return self.check_db_state()

    def check_db_state(self):
        # Copy over our db to a named temp file so we can open it
        with NamedTemporaryFile() as db_ntf:
            with self.target.open('rb') as db_h:
                db_ntf.write(db_h.read())
            con = sqlite3.connect(db_ntf.name)
            cursor = con.cursor()
            # Check for the expected tables
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = {t[0] for t in cursor.fetchall()}
            expected_tables = {
                'multiclass',
                'placements',
                'ranks',
                'taxa'
            }
            # We can tolerate extra tables, but at least these should exist
            if not tables.issuperset(expected_tables):
                log.warn("Placement DB present but missing tables {}".format(
                    (expected_tables - tables)
                ))
                if self.removed_malformed:
                    log.warn("Removing malformed Placement DB")
                    self.target.remove()
                return False
            # Implicit else we have tables....
            if self.expected_tax_ids:
                cursor.execute("SELECT DISTINCT tax_id from taxa;")
                db_taxa = {t[0] for t in cursor.fetchall()}
                if not db_taxa.issuperset(self.expected_tax_ids):
                    log.warn("DB exists, but missing tax IDs {}".format(
                        self.expected_tax_ids - db_taxa
                    ))
                    if self.removed_malformed:
                        log.warn("Removing malformed Placement DB")
                        self.target.remove()
                    return False
            if self.expected_ranks:
                cursor.execute("SELECT DISTINCT rank from ranks;")
                db_ranks = {t[0] for t in cursor.fetchall()}
                if not db_ranks.issuperset(self.expected_ranks):
                    log.warn("DB exists, but missing ranks {}".format(
                        self.expected_ranks - db_ranks
                    ))
                    if self.removed_malformed:
                        log.warn("Removing malformed Placement DB")
                        self.target.remove()
                    return False
            # Implicit else we passed all our tests...
            return True


class PlacementDB_Classified_ContainerTargetInfo(sl.ContainerTargetInfo):
    def __init__(
            self,
            task,
            path,
            is_tmp=False,
            client=None,
            guppy_command=None,
            ):
                super(PlacementDB_Classified_ContainerTargetInfo, self).__init__(
                    task=task,
                    path=path,
                    format=luigi.format.Nop,
                    is_tmp=is_tmp,
                    client=client
                    )
                self.target.basic_exists = self.target.exists
                self.target.exists = self.full_exists
                self.guppy_command = guppy_command
                self.guppy_namespace = self.guppy_classify_parser().parse_known_args(
                    shlex.split(guppy_command))[0]

    def guppy_classify_parser(self):
        parser = argparse.ArgumentParser()
        parser.add_argument(
            '--seed',
            default=1,
        )
        parser.add_argument(
            '--classifier',
            default='pplacer'
        )
        parser.add_argument(
            '--cutoff',
            default=0.9
        )
        parser.add_argument(
            '--bayes-cutoff',
            default=1.00
        )
        parser.add_argument(
            '--multiclass-min',
            default=0.2
        )
        parser.add_argument(
            '--bootstrap-cutoff',
            default=0.8
        )
        parser.add_argument(
            '--bootstrap-extension-cutoff',
            default=0.4
        )
        parser.add_argument(
            '--pp',
            action='store_true'
        )
        parser.add_argument(
            '--word-length',
            default=8
        )
        parser.add_argument(
            '--nbc-rank',
            default='genus'
        )
        parser.add_argument(
            '--n-boot',
            default=100,
        )
        parser.add_argument(
            '--no-pre-mask',
            action='store_true'
        )
        parser.add_argument(
            '--nbc-as-rdp',
            action='store_true'
        )
        parser.add_argument(
            '--no-random-tie-break',
            action='store_true'
        )
        return parser

    def full_exists(self):
        if not self.target.basic_exists():
            return False
        # Implicit else, the file exists. Check our db for completeness
        return self.check_db_state()

    def check_db_state(self):
        # Copy over our db to a named temp file so we can open it
        with NamedTemporaryFile() as db_ntf:
            with self.target.open('rb') as db_h:
                db_ntf.write(db_h.read())
            con = sqlite3.connect(db_ntf.name)
            cursor = con.cursor()
            # Get all the classification runs stored in this database
            cursor.execute("SELECT params from runs;")
            completed_guppy_namespace = [
                self.guppy_classify_parser().parse_known_args(
                    shlex.split(p[0]))[0]
                for p in cursor.fetchall()]
            if not self.guppy_command:
                # No classifier specified. If we have any runs, assume we are OK
                if len(completed_guppy_namespace) > 0:
                    return True
            elif self.guppy_namespace not in completed_guppy_namespace:
                return False
            # Implicit else all tests passed
            return True


class PlacementDB_MCC_ContainerTargetInfo(sl.ContainerTargetInfo):
    def __init__(
            self,
            task,
            path,
            is_tmp=False,
            client=None,
            expected_seq_ids=None):
        super(PlacementDB_MCC_ContainerTargetInfo, self).__init__(
            task=task,
            path=path,
            format=luigi.format.Nop,
            is_tmp=is_tmp,
            client=client
            )
        self.target.basic_exists = self.target.exists
        self.target.exists = self.full_exists
        self.expected_seq_ids = expected_seq_ids

    def full_exists(self):
        if not self.target.basic_exists():
            return False
        # Implicit else, the file exists. Check our db for completeness
        return self.check_db_state()

    def check_db_state(self):
        # Copy over our db to a named temp file so we can open it
        with NamedTemporaryFile() as db_ntf:
            with self.target.open('rb') as db_h:
                db_ntf.write(db_h.read())
            con = sqlite3.connect(db_ntf.name)
            cursor = con.cursor()
            # Get the tables
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = {t[0] for t in cursor.fetchall()}
            if 'multiclass_concat' not in tables:
                return False
            if self.expected_seq_ids:
                # We have some seq_ids. Make sure we have something in mcc for each
                cursor.execute("SELECT DISTINCT name FROM multiclass_concat;")
                db_seq_ids = {r[0] for r in cursor.fetchall()}
                if not db_seq_ids.issuperset(self.expected_seq_ids):
                    log.warn("Multiclass concat is missing seq ids {}".format(
                        self.expected_seq_ids - db_seq_ids
                    ))
                    return False

            # Implicit else...
            return True


class PlacementDB_SI_ContainerTargetInfo(sl.ContainerTargetInfo):
    def __init__(
            self,
            task,
            path,
            is_tmp=False,
            client=None,
            expected_seq_ids=None):
        super(PlacementDB_SI_ContainerTargetInfo, self).__init__(
            task=task,
            path=path,
            format=luigi.format.Nop,
            is_tmp=is_tmp,
            client=client
            )
        self.target.basic_exists = self.target.exists
        self.target.exists = self.full_exists

        self.expected_seq_ids = expected_seq_ids

    def full_exists(self):
        if not self.target.basic_exists():
            return False
        # Implicit else, the file exists. Check our db for completeness
        return self.check_db_state()

    def check_db_state(self):
        # Copy over our db to a named temp file so we can open it
        with NamedTemporaryFile() as db_ntf:
            with self.target.open('rb') as db_h:
                db_ntf.write(db_h.read())
            con = sqlite3.connect(db_ntf.name)
            cursor = con.cursor()
            # Get the tables
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = {t[0] for t in cursor.fetchall()}
            if 'seq_info' not in tables:
                return False
            if self.expected_seq_ids:
                # If we recieved what the expected seq ids are
                cursor.execute("SELECT name FROM seq_info;")
                db_seq_ids = {r[0] for r in cursor.fetchall()}
                if not db_seq_ids.issuperset(self.expected_seq_ids):
                    log.warn("Seq Info exists, but is missing {}".format(
                        self.expected_seq_ids - db_seq_ids
                    ))
                    # Missing some entries. Drop the table and rebuild it
                    cursor.executescript("DROP TABLE seq_info;")
                    con.commit()
                    # Now have to copy the db back to the source
                    with self.target.open('wb') as db_h:
                        db_h.write(db_ntf.read())
                    return False

            # Implicit else...
            return True
