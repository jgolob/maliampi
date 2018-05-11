import sciluigi as sl
import luigi
from .NCBI_NT_Repository import NCBI_NT_Repository

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
# Check if every entry has a feature table.


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
