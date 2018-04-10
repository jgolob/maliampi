import luigi
import sciluigi as sl
import os
from string import Template
from Bio import AlignIO
import shutil


# Tasks
class LoadFastaSeqs(sl.ExternalTask):
    fasta_seq_path = sl.Parameter()

    def out_seqs(self):
        return sl.TargetInfo(self, self.fasta_seq_path)


class LoadFile(sl.ExternalTask):
    path = sl.Parameter()

    def out(self):
        return sl.TargetInfo(self, self.path)


class SearchRepoForMatches(sl.ContainerTask):
    # A Task that uses vsearch to find matches for experimental sequences in a repo of sequences

    # Define the container (in docker-style repo format) to complete this task
    container = 'docker://golob/vsearch'

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
        # Get our host paths for inputs and outputs
        input_paths = {
            'exp_seqs': self.in_exp_seqs().path,
            'repo_seqs': self.in_repo_seqs().path,
        }

        print(self.out_matched_repo_seqs().path)
        output_paths = {
            'matched_repo': self.out_matched_repo_seqs().path,
            'unmatched_exp': self.out_unmatched_exp_seqs().path,
            'uc': self.out_matches_uc().path,
        }

        self.ex(
            command='vsearch ' +
                    ' --threads=%s' % self.containerinfo.vcpu +
                    ' --usearch_global $exp_seqs' +
                    ' --db=$repo_seqs' +
                    ' --id=%s' % self.min_id +
                    ' --strand both' +
                    ' --uc=$uc --uc_allhits' +
                    ' --notmatched=$unmatched_exp' +
                    ' --dbmatched=$matched_repo' +
                    ' --maxaccepts=%s' % self.maxaccepts,
            input_paths=input_paths,
            output_paths=output_paths,
            )


class FillLonely(sl.Task):
    pass


class CMAlignSeqs(sl.ContainerTask):
    # A Task that uses CMAlign to make an alignment
    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/infernal'

    in_seqs = None  # Seqs to align

    # Where to put the alignment (both sto and fasta created) include prefix but not file extension

    alignment_sto_fn = sl.Parameter()

    # FULL path where to store the alignment scores. Include extension
    alignment_score_fn = sl.Parameter()

    # cmalign parameters
    # max memory to use in MB.
    cmalign_mxsize = sl.Parameter(default=8196)

    def out_alignscores(self):
        return sl.TargetInfo(self, self.alignment_score_fn)

    def out_align_sto(self):
        return sl.TargetInfo(self, self.alignment_sto_fn)

    def run(self):

        # Get our host paths for inputs and outputs
        input_paths = {
            'in_seqs': self.in_seqs().path,
        }

        output_paths = {
            'align_sto': self.out_align_sto().path,
            'alignscores': self.out_alignscores().path,
        }

        self.ex(
            command='cmalign' +
            ' --cpu %s --noprob --dnaout --mxsize %d' % (self.containerinfo.vcpu, self.cmalign_mxsize) +
            ' --sfile $alignscores -o $align_sto ' +
            ' /cmalign/data/SSU_rRNA_bacteria.cm $in_seqs',
            input_paths=input_paths,
            output_paths=output_paths
        )


class AlignmentStoToFasta(sl.Task):
    # Alignment in stockholm format
    in_align_sto = None

    align_fasta_fn = sl.Parameter()

    def out_align_fasta(self):
        return sl.TargetInfo(self, self.align_fasta_fn)

    def run(self):
        # Use biopython to convert from stockholm to fasta output
        AlignIO.write(AlignIO.read(self.in_align_sto().open(), 'stockholm'), self.out_align_fasta().path, 'fasta')


class RAxMLTree(sl.ContainerTask):
    # A task that uses RAxML to generate a tree from an alignment

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/raxml'

    # Input of an alignment in FASTA format
    in_align_fasta = None

    # Parameter: Path + filename where the resultant tree should go
    tree_path = sl.Parameter()
    tree_stats_path = sl.Parameter()
    # DIRECTORY where the intermediate RAxML files should go
    raxml_working_dir = sl.Parameter()
    
    # Parameters for RAxML

    raxml_model = sl.Parameter(default='GTRGAMMA')
    raxml_parsimony_seed = sl.Parameter(default=12345)

    def out_tree(self):
        return sl.TargetInfo(self, self.tree_path)

    def out_tree_stats(self):
        return sl.TargetInfo(self, self.tree_stats_path)

    def run(self):
        # Lots of filesystem throat-clearing
        name = os.path.basename(os.path.splitext(self.tree_path)[0])
        if not os.path.exists(self.raxml_working_dir):
            os.makedirs(self.raxml_working_dir)
        raxml_working_dir = os.path.abspath(self.raxml_working_dir)

        # Get our host paths for inputs and outputs
        # To be mapped into the container as appropriate
        input_paths = {
            'in_align_fasta': self.in_align_fasta().path,
        }

        output_paths = {
            'out_tree': self.out_tree().path,
            'out_tree_stats': self.out_tree_stats().path,
            'raxml_working_dir': self.raxml_working_dir,
        }

        self.ex(
            command='raxml'
                    ' -n %s' % name +  # Prefix/name to use for the output files
                    ' -m %s' % self.raxml_model +  # Model to use
                    ' -s $in_align_fasta' +  # Path to input alignment
                    ' -p %d' % self.raxml_parsimony_seed +
                    ' -T %d' % max([int(self.containerinfo.vcpu), 2]) +  # threads (min 2)
                    ' -w $raxml_working_dir',  # working directory
            input_paths=input_paths,
            output_paths=output_paths,
            inputs_mode='rw',
        )

        # Move the resultant tree to our specified path
        shutil.copy(os.path.join(raxml_working_dir, "RAxML_bestTree.%s" % name), self.out_tree().path)
        # And stats
        shutil.copy(os.path.join(raxml_working_dir, "RAxML_info.%s" % name), self.out_tree_stats().path)
