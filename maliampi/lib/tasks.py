import luigi
import sciluigi as sl
import os
from string import Template
from Bio import AlignIO
import shutil
from .targets import NCBI_Repo_Entries_TargetInfo,  NCBI_Repo_Filled_TargetInfo
from lib.ExtractGenbank import ExtractGenbank
import csv
import uuid
import json
from collections import defaultdict
import logging
from datetime import datetime
import pytz

log = logging.getLogger('sciluigi-interface')


# Tasks
class LoadFastaSeqs(sl.ExternalTask):
    fasta_seq_path = sl.Parameter()

    def out_seqs(self):
        return sl.ContainerTargetInfo(self, self.fasta_seq_path)


class LoadFile(sl.ExternalTask):
    path = sl.Parameter()

    def out_file(self):
        return sl.ContainerTargetInfo(self, self.path)


class SearchRepoForMatches(sl.ContainerTask):
    # A Task that uses vsearch to find matches for experimental sequences in a repo of sequences

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/vsearch:2.7.1_bcw_0.2.0'

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
        return sl.ContainerTargetInfo(self, self.matches_uc_path)

    def out_unmatched_exp_seqs(self):
        return sl.ContainerTargetInfo(self, self.unmatched_exp_seqs_path)

    def out_matched_repo_seqs(self):
        return sl.ContainerTargetInfo(self, self.matched_repo_seqs_path)

    def run(self):
        # Get our host paths for inputs and outputs
        input_targets = {
            'exp_seqs': self.in_exp_seqs(),
            'repo_seqs': self.in_repo_seqs(),
        }
        output_targets = {
            'matched_repo': self.out_matched_repo_seqs(),
            'unmatched_exp': self.out_unmatched_exp_seqs(),
            'uc': self.out_matches_uc(),
        }

        self.ex(
            command='vsearch ' +
                    ' --threads=%s' % self.containerinfo.vcpu +
                    ' --usearch_global $exp_seqs' +
                    ' --db $repo_seqs' +
                    ' --id=%s' % self.min_id +
                    ' --strand both' +
                    ' --uc=$uc --uc_allhits' +
                    ' --notmatched=$unmatched_exp' +
                    ' --dbmatched=$matched_repo' +
                    ' --maxaccepts=%s' % self.maxaccepts,
            input_targets=input_targets,
            output_targets=output_targets,
            )


class FillLonely(sl.Task):
    pass


class CMAlignSeqs(sl.ContainerTask):
    # A Task that uses CMAlign to make an alignment
    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/infernal:1.1.2_bcw_0.2.0'

    in_seqs = None  # Seqs to align

    # Where to put the alignment (both sto and fasta created) include prefix but not file extension

    alignment_sto_fn = sl.Parameter()

    # FULL path where to store the alignment scores. Include extension
    alignment_score_fn = sl.Parameter()

    # cmalign parameters
    # max memory to use in MB.
    cmalign_mxsize = sl.Parameter(default=8196)

    def out_alignscores(self):
        return sl.ContainerTargetInfo(self, self.alignment_score_fn)

    def out_align_sto(self):
        return sl.ContainerTargetInfo(self, self.alignment_sto_fn)

    def run(self):

        # Get our host paths for inputs and outputs
        input_targets = {
            'in_seqs': self.in_seqs(),
        }

        output_targets = {
            'align_sto': self.out_align_sto(),
            'alignscores': self.out_alignscores(),
        }

        self.ex(
            command='cmalign' +
            ' --cpu %s --noprob --dnaout --mxsize %d' % (self.containerinfo.vcpu, self.cmalign_mxsize) +
            ' --sfile $alignscores -o $align_sto ' +
            ' /cmalign/data/SSU_rRNA_bacteria.cm $in_seqs',
            input_targets=input_targets,
            output_targets=output_targets
        )


class AlignmentStoToFasta(sl.Task):
    # Alignment in stockholm format
    in_align_sto = None

    align_fasta_fn = sl.Parameter()

    def out_align_fasta(self):
        return sl.ContainerTargetInfo(self, self.align_fasta_fn)

    def run(self):
        # Use biopython to convert from stockholm to fasta output
        with self.out_align_fasta().open('w') as out_h:
            AlignIO.write(
                        AlignIO.read(self.in_align_sto().open('r'), 'stockholm'),
                        out_h,
                        'fasta'
            )


class RAxMLTree(sl.ContainerTask):
    # A task that uses RAxML to generate a tree from an alignment

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/raxml:8.2.11_bcw_0.2.0'

    # Input of an alignment in FASTA format
    in_align_fasta = None

    # Parameter: Path + filename where the resultant tree should go
    tree_path = sl.Parameter()
    tree_stats_path = sl.Parameter()
    # DIRECTORY where the intermediate RAxML files should go (container fs space)
    raxml_working_dir = sl.Parameter(default='/scratch')

    # Parameters for RAxML

    raxml_model = sl.Parameter(default='GTRGAMMA')
    raxml_parsimony_seed = sl.Parameter(default=12345)

    def out_tree(self):
        return sl.ContainerTargetInfo(self, self.tree_path)

    def out_tree_stats(self):
        return sl.ContainerTargetInfo(self, self.tree_stats_path)

    def run(self):
        # Lots of filesystem throat-clearing
        name = os.path.basename(os.path.splitext(self.tree_path)[0])

        # Get our host paths for inputs and outputs
        # To be mapped into the container as appropriate
        input_targets = {
            'in_align_fasta': self.in_align_fasta(),
        }

        output_targets = {
            'out_tree': self.out_tree(),
            'out_tree_stats': self.out_tree_stats(),
        }

        self.ex(
            command='mkdir -p $raxml_working_dir && raxml'
                    ' -n %s' % name +  # Prefix/name to use for the output files
                    ' -m %s' % self.raxml_model +  # Model to use
                    ' -s $in_align_fasta' +  # Path to input alignment
                    ' -p %d' % self.raxml_parsimony_seed +
                    ' -T %d' % max([int(self.containerinfo.vcpu), 2]) +  # threads (min 2)
                    ' -w $raxml_working_dir' +  # working directory
                    ' && cp $raxml_working_dir/RAxML_bestTree.%s $out_tree' % name +
                    ' && cp $raxml_working_dir/RAxML_info.%s $out_tree_stats' % name,
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={'raxml_working_dir': self.raxml_working_dir},
            inputs_mode='rw',
        )


class NT_AccessionsForQuery(sl.ContainerTask):
    # A task that takes a query and returns all of the NCBI NT accessions
    # matching
    container = 'golob/medirect:0.9.0_BCW_0.2.0a'

    query = sl.Parameter()
    email = sl.Parameter()
    ncbi_concurrent_connections = sl.Parameter(default=3)
    retry_max = sl.Parameter(default=1)
    retry_delay = sl.Parameter(default=60000)

    accessions_fn = sl.Parameter()

    def out_accessions(self):
        return sl.ContainerTargetInfo(self, self.accessions_fn)

    def run(self):
        if (self.out_accessions().scheme == 'file'):
            os.makedirs(os.path.dirname(self.out_accessions().path))
        output_targets = {
            'accessions': self.out_accessions()
        }

        self.ex(
            command='ncbi_get_nt_accessions_for_query' +
                    ' --email %s' % self.email +
                    ' --query "%s"' % self.query +
                    ' --ncbi_concurrent_connections %d ' % self.ncbi_concurrent_connections +
                    ' --retry_max %d' % self.retry_max +
                    ' --retry_delay %d' % self.retry_delay +
                    ' --out $accessions',
            output_targets=output_targets
        )


class NT_Repo_Update_Accessions(sl.Task):
    # Takes a file with accesions. Inserts new entries into the repo
    in_accessions = None
    in_repo_url = None

    extra_values = sl.Parameter(default={})

    def get_accessions(self):
        if not self.in_accessions().target.exists():
            return -1
        else:
            return {
                r[0] for r
                in csv.reader(self.in_accessions().open())
            }

    def out_repo(self):
        return NCBI_Repo_Entries_TargetInfo(
            self,
            self.in_repo_url().open().read().strip(),
            self.get_accessions()
            )

    def run(self):
        self.out_repo().target.repo.add_new_versions(
            self.get_accessions(),
            extra_values=json.loads(self.extra_values)
            )


class NT_Repo_Fill(sl.ContainerTask):
    """Finds skeletal entries in our repo, and retrieves the entry from genbank
    """
    container = 'golob/medirect:0.9.0_BCW_0.2.0a'
    # Takes a file with accesions. Inserts new entries into the repo
    in_repo = None

    constraints = sl.Parameter(default={})
    debug = sl.Parameter(default=False)
    ncbi_threads_after_hours = sl.Parameter(default=10)
    ncbi_threads_peak_hours = sl.Parameter(default=3)
    working_dir = sl.Parameter()
    email = sl.Parameter()

    chunk_size = sl.Parameter(default=100)

    def out_repo(self):
        return NCBI_Repo_Filled_TargetInfo(
            self,
            self.in_repo().path,
            )

    def chunks(self, l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def work_on_chunk(self, chunk_versions, gb_needed, raw_gb):
        ncbi_threads = 3
        ncbi_tz = pytz.timezone('US/Eastern')
        with gb_needed.open('w') as gb_needed_h:
                gb_needed_w = csv.writer(gb_needed_h)
                gb_needed_w.writerow(['id'])
                for v in chunk_versions:
                    gb_needed_w.writerow([v])

        # Figure out if we are after hours for NCBI
        ncbi_time = datetime.now(ncbi_tz)
        # Is it a weekend?
        if ncbi_time.weekday() >= 5:
            ncbi_threads = self.ncbi_threads_after_hours
        elif ((ncbi_time.hour <= 18) and (ncbi_time.hour > 6)):
            # It's a weekday, and between 6 am - 6pm, not after hours
            ncbi_threads = self.ncbi_threads_peak_hours
        else:
            # After hours on a weekday
            ncbi_threads = self.ncbi_threads_after_hours

        self.ex(
            command="mefetch -id $gb_needed" +
                    " -db nucleotide -format gbwithparts -mode text -csv -retmax 1" +
                    " -email %s" % self.email +
                    " -proc %d" % int(ncbi_threads) +
                    " -out $raw_gb",
            input_targets={'gb_needed': gb_needed},
            output_targets={'raw_gb': raw_gb}
        )

        records = ExtractGenbank(raw_gb.open('r'))

        for version, record in records.get_records().items():
            logging.info("Updating document for {}".format(version))
            self.out_repo().target.repo.add_entry(record)

    def run(self):
        # Determine when this is running to see how hard we can hit NCBI
        versions_need_fill = self.out_repo().target.repo.versions_needing_data(
            self.constraints
        )
        versions_need_fill = list(versions_need_fill)
        if self.debug:
            versions_need_fill = versions_need_fill[0:2]
        random_dir = str(uuid.uuid4())
        gb_needed = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.working_dir,
                'gb',
                random_dir,
                'needed.csv'
            ))
        raw_gb = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.working_dir,
                'gb',
                random_dir,
                'raw_gb.csv'
            ))
        for chunk_i, chunk_versions in enumerate(self.chunks(
                versions_need_fill,
                int(self.chunk_size))
                ):
            log.info("Updating chunk {} of {}".format(
                chunk_i+1,
                int(len(versions_need_fill) / self.chunk_size))
            )
            self.work_on_chunk(chunk_versions, gb_needed, raw_gb)
