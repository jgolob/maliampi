import luigi
import sciluigi as sl
from lib.containertask import ContainerTask, ContainerInfo
import os
from string import Template


# Tasks
class LoadFastaSeqs(sl.ExternalTask):
    fasta_seq_path = sl.Parameter()

    def out_seqs(self):
        return sl.TargetInfo(self, self.fasta_seq_path)


class ListVsearchOptions(ContainerTask):
    # A Task that uses vsearch to find matches for experimental sequences in a repo of sequences

    # Define the container (in docker-style repo format) to complete this task
    __container__ = 'golob/vsearch'

    def out_hello(self):
        return sl.TargetInfo(self, "../working/hello.txt")

    def run(self):

        self.ex(
                'vsearch' +
                ' -h'
            )


class SearchRepoForMatches(ContainerTask):
    # A Task that uses vsearch to find matches for experimental sequences in a repo of sequences

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/vsearch'
    num_cpu = 1

    in_exp_seqs = None  # Experimental seqs for this task
    in_repo_seqs = None  # Repository seqs

    # Where to put the matches in the repo, in UC format
    matches_uc_path = sl.Parameter()
    # Where to put the unmatched exp seqs
    unmatched_exp_seqs_path = sl.Parameter()
    # Where to put the repo seqs with at least one match in the exp
    matched_repo_seqs_path = sl.Parameter()

    # vsearch parameters
    min_id = sl.Parameter(default=1.0)
    maxaccepts = sl.Parameter(default=1)  # by default, stop searching after the first

    def out_matches_uc(self):
        return sl.TargetInfo(self, self.matches_uc_path)

    def out_unmatched_exp_seqs(self):
        return sl.TargetInfo(self, self.unmatched_exp_seqs_path)

    def out_matched_repo_seqs(self):
        return sl.TargetInfo(self, self.matched_repo_seqs_path)

    def run(self):
        input_paths = {
            'exp_seqs': self.in_exp_seqs().path,
            'repo_seqs': self.in_repo_seqs().path,
        }

        output_paths = {
            'matched_repo': self.matched_repo_seqs_path,
            'unmatched_exp': self.unmatched_exp_seqs_path,
            'uc': self.matches_uc_path,
        }

        self.ex(
            command='vsearch '+
                  ' --threads=%s' % self.num_cpu+
                  ' --usearch_global $exp_seqs'+
                  ' --db=$repo_seqs'+
                  ' --id=%s' % self.min_id+
                  ' --strand both'+
                  ' --uc=$uc --uc_allhits'+
                  ' --notmatched=$unmatched_exp'+
                  ' --dbmatched=$matched_repo'+
                  ' --maxaccepts=%s' % self.maxaccepts,
            input_paths=input_paths,
            output_paths=output_paths,
            )
