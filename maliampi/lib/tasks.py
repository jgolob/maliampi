import luigi
import sciluigi as sl
from lib.containertask import ContainerTask, ContainerInfo
import os


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
    __container__ = 'golob/vsearch'
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
        input_paths = [
            self.in_exp_seqs().path,
            self.in_repo_seqs().path,
        ]
        output_paths = [
            self.matched_repo_seqs_path,
            self.unmatched_exp_seqs_path,
            self.matches_uc_path,
        ]
        input_common_prefix = os.path.commonprefix(input_paths)
        output_common_prefix = os.path.commonprefix(output_paths)

        command = 'vsearch' + \
                ' --threads={}'.format(self.num_cpu) + \
                ' --usearch_global {}'.format(
                        os.path.relpath(
                            self.in_exp_seqs().path,
                            input_common_prefix)) + \
                ' --db={}'.format(os.path.relpath(
                    self.in_repo_seqs().path,
                    input_common_prefix)) + \
                ' --id={}'.format(self.min_id) + \
                ' --strand both' + \
                ' --uc={} --uc_allhits'.format(os.path.relpath(
                    self.out_matches_uc().path,
                    output_common_prefix)) + \
                ' --notmatched={}'.format(os.path.relpath(
                    self.out_unmatched_exp_seqs().path,
                    output_common_prefix)) + \
                ' --dbmatched={}'.format(os.path.relpath(
                        self.out_matched_repo_seqs().path,
                        output_common_prefix)) + \
                ' --maxaccepts={}'.format(self.maxaccepts)

        print(command)
        print(input_common_prefix)
        print(output_common_prefix)

        self.ex(command=
                'vsearch' +
                ' --threads={}'.format(self.num_cpu) +
                ' --usearch_global {}'.format(
                    os.path.relpath(self.in_exp_seqs().path, input_common_prefix)) +
                ' --db={}'.format(os.path.relpath(self.in_repo_seqs().path, input_common_prefix)) +
                ' --id={}'.format(self.min_id) +
                ' --strand both' +
                ' --uc={} --uc_allhits'.format(os.path.relpath(
                    self.out_matches_uc().path,
                    output_common_prefix)) +
                ' --notmatched={}'.format(os.path.relpath(
                    self.out_unmatched_exp_seqs().path,
                    output_common_prefix)) +
                ' --dbmatched={}'.format(os.path.relpath(
                        self.out_matched_repo_seqs().path,
                        output_common_prefix)) +
                ' --maxaccepts={}'.format(self.maxaccepts),
                mounts={
                    input_common_prefix: {'bind': '/mnt/input', 'mode': 'ro'},
                    output_common_prefix: {'bind': '/mnt/output', 'mode': 'rw'}
                }

            )
