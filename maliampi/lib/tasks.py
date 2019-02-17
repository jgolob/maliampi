import luigi
import sciluigi as sl
import os
from string import Template
from Bio import AlignIO, SeqIO
import shutil
from .targets import NCBI_Repo_Entries_TargetInfo,  NCBI_Repo_Filled_TargetInfo
from .targets import PlacementDB_Prepped_ContainerTargetInfo, RefpkgTGZ_ContainerTargetInfo
from .targets import PlacementDB_Classified_ContainerTargetInfo, NCBI_Repo_Peptides_TargetInfo
from .targets import PlacementDB_MCC_ContainerTargetInfo, PlacementDB_SI_ContainerTargetInfo
from lib.ExtractGenbank import ExtractGenbank
import csv
import uuid
import json
from collections import defaultdict
import logging
from datetime import datetime
import pytz
import tarfile
import io
import gzip
import json
import gc
import tempfile
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool

log = logging.getLogger('sciluigi-interface')


# Tasks
class LoadFastaSeqs(sl.ExternalTask):
    fasta_seq_path = sl.Parameter()

    def out_seqs(self):
        return sl.ContainerTargetInfo(self, self.fasta_seq_path)


class LoadFile(sl.ExternalTask):
    path = sl.Parameter()
    file_format = sl.Parameter(default=None)

    def out_file(self):
        if self.file_format == 'gzip':
            file_format = luigi.format.Gzip
        else:
            file_format = None

        return sl.ContainerTargetInfo(self, self.path, format=file_format)


class LoadRefpkgTGZ(sl.ExternalTask):
    path = sl.Parameter()

    def out_refpkg_tgz(self):
        return RefpkgTGZ_ContainerTargetInfo(self, self.path)


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


class VerifyRepo(sl.ContainerTask):
    # Verify repo sequences are somewhat like expected sequences

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/vsearch:2.7.1_bcw_0.2.0'

    in_repo_seqs = None  # Repo seqs
    in_expected_seqs = None  # Expected seqs

    # Where to put the matches in the repo, in UC format
    uc_fn = sl.Parameter()

    # Where to put the unmatched exp seqs
    unverified_seqs_fn = sl.Parameter()
    # Where to put the repo seqs with at least one match in the exp
    verified_seqs_fn = sl.Parameter()

    # vsearch parameters
    min_id = sl.Parameter(default=0.7)
    maxaccepts = sl.Parameter(default=1)  # by default, stop searching after the first
    query_cov = sl.Parameter(default=0.70)
    iddef = sl.Parameter(default=2)

    def out_repo_uc(self):
        return sl.ContainerTargetInfo(self, self.uc_fn)

    def out_verified_seqs(self):
        return sl.ContainerTargetInfo(self, self.verified_seqs_fn)

    def out_unverified_seqs(self):
        return sl.ContainerTargetInfo(self, self.unverified_seqs_fn)

    def run(self):
        # Get our host paths for inputs and outputs
        input_targets = {
            'repo_seqs': self.in_repo_seqs(),
            'expected_seqs': self.in_expected_seqs(),
        }
        output_targets = {
            'verified': self.out_verified_seqs(),
            'unverified': self.out_unverified_seqs(),
            'uc': self.out_repo_uc(),
        }

        self.ex(
            command='vsearch ' +
                    ' --threads=%s' % self.containerinfo.vcpu +
                    ' --usearch_global $repo_seqs' +
                    ' --db $expected_seqs' +
                    ' --id=$min_id' +
                    ' --iddef $iddef ' +
                    ' --query_cov $query_cov ' +
                    ' --strand both' +
                    ' --uc=$uc ' +
                    ' --notmatched=$unverified' +
                    ' --matched=$verified ' +
                    ' --maxaccepts=%s' % self.maxaccepts,
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'min_id': self.min_id,
                'iddef': self.iddef,
                'query_cov': self.query_cov
                }
            )


class CombineRepoMatches(sl.Task):
    # Given input LISTS of recruited sequences, combine to one master list
    # AND filter seq info to match
    in_seqs = None
    in_seq_info = None

    # Where to put the combined seqs
    seqs_fn = sl.Parameter()
    # Where to put the combined seq info
    seq_info_fn = sl.Parameter()

    def out_seqs(self):
        return sl.ContainerTargetInfo(self, self.seqs_fn)

    def out_seq_info(self):
        return sl.ContainerTargetInfo(self, self.seq_info_fn)

    def run(self):
        seqs_included = set()
        seq_ids_included = set()
        with self.out_seqs().open('w') as out_seq_h:
            for seq_i, seq_T in enumerate(self.in_seqs):
                if seq_i == 0:  # We are on the first block. Include everything
                    with seq_T().open('r') as seq_h:
                        for sr in SeqIO.parse(seq_h, 'fasta'):
                            seqs_included.add(sr.seq)
                            seq_ids_included.add(sr.id)
                            SeqIO.write(sr, out_seq_h, 'fasta')
                else:  # AFTER the first block
                    with seq_T().open('r') as seq_h:
                        for sr in SeqIO.parse(seq_h, 'fasta'):
                            if sr.id in seq_ids_included:
                                continue
                            # Implicit else
                            if sr.seq not in seqs_included:
                                continue
                            # implicit else
                            seqs_included.add(sr.seq)
                            seq_ids_included.add(sr.id)
                            SeqIO.write(sr, out_seq_h, 'fasta')

        log.info("{} combined sequences".format(len(seq_ids_included)))
        # Now work on the seq_info
        # First figure out the full set of columns for both seq_info tables
        seq_info_headers = [
            csv.DictReader(T().open('r')).fieldnames
            for T in self.in_seq_info
        ]
        shared_headers = []
        for h_i, h in enumerate(seq_info_headers):
            if h_i == 0:
                shared_headers = h
            else:  # Beyond the first
                for c in h:
                    if c not in shared_headers:
                        shared_headers.append(c)

        with self.out_seq_info().open('w') as out_seq_info_h:
            si_writer = csv.DictWriter(
                out_seq_info_h,
                fieldnames=shared_headers
                )
            si_writer.writeheader()
            for si_T in self.in_seq_info:
                with si_T().open('r') as si_h:
                    si_R = csv.DictReader(si_h)
                    for row in si_R:
                        if row.get('seqname') in seq_ids_included:
                            si_writer.writerow({
                                c: row.get(c, "")
                                for c in shared_headers
                            })
                            seq_ids_included.remove(row.get('seqname'))


class FillLonely(sl.Task):
    pass


class CMSearchVerify(sl.ContainerTask):
    # A Task that uses cmsearch to check if seqs are valid
    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/infernal:1.1.2_bcw_0.2.0'

    in_seqs = None  # Seqs to verify

    # Where to put the table of search results

    results_fn = sl.Parameter()

    def out_results(self):
        return sl.ContainerTargetInfo(
            self,
            self.results_fn
            )

    def run(self):

        # Get our host paths for inputs and outputs
        input_targets = {
            'in_seqs': self.in_seqs(),
        }

        output_targets = {
            'results': self.out_results(),
        }

        self.ex(
            command=(
                'cmsearch '
                '--cpu $vcpu '
                '--noali '
                '--tblout $results '
                '/cmalign/data/SSU_rRNA_bacteria.cm '
                '$in_seqs'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vcpu': self.containerinfo.vcpu,
            }
        )


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
        working_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4())
        )
        # Get our host paths for inputs and outputs
        input_targets = {
            'in_seqs': self.in_seqs(),
        }

        output_targets = {
            'align_sto': self.out_align_sto(),
            'alignscores': self.out_alignscores(),
        }

        self.ex(
            command=(
                'mkdir -p $working_dir && cd $working_dir &&'
                ' cmalign'
                ' --cpu $vcpu --noprob --dnaout --mxsize $mxsize' 
                ' --sfile $alignscores -o $align_sto '
                ' /cmalign/data/SSU_rRNA_bacteria.cm $in_seqs'
                ' && rm -r $working_dir'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'vcpu': self.containerinfo.vcpu,
                'mxsize': self.cmalign_mxsize,
                'working_dir': working_dir
            }
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


class CleanupTreeInfo(sl.Task):
    #  Get rid of extraneous info from tree info file from RAxML
    in_tree_info = None
    tree_info_path = sl.Parameter()

    def out_tree_info(self):
        return sl.ContainerTargetInfo(
            self,
            self.tree_info_path
        )

    def run(self):
        with self.out_tree_info().open('w') as out_h:
            with self.in_tree_info().open('r') as in_h:
                past_cruft = False
                for l in in_h:
                    if "This is RAxML version" == l[0:21]:
                        past_cruft = True
                    if past_cruft:
                        out_h.write(l)


class ConfirmSeqInfoTaxonomy(sl.ContainerTask):
    #  tax ids can change / be retired over time.
    #  Confirm all tax ids in our seq info file are current
    #  Correct as needed
    container = 'golob/seqinfo_taxonomy_sync:0.2.1__bcw.0.3.0'

    in_tax_db = None
    in_seq_info = None

    # Entrez expects an email
    email = sl.Parameter()
    # Where to put the corrected seqinfo
    confirmed_seqinfo_path = sl.Parameter()

    def out_seq_info(self):
        return sl.ContainerTargetInfo(
            self,
            self.confirmed_seqinfo_path
        )
    
    def run(self):
        input_targets = {
            'seqinfo': self.in_seq_info(),
            'taxonomy_db': self.in_tax_db(),
        }
        output_targets = {
            'outseqinfo': self.out_seq_info()
        }
        extra_params = {
            'email': self.email
        }
        self.ex(
            command=(
                'seqinfo_taxonomy_sync.py '
                '$seqinfo '
                '$outseqinfo '
                ' --db $taxonomy_db '
                ' --email $email'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params=extra_params
        )


class RAxMLTree(sl.ContainerTask):
    # A task that uses RAxML to generate a tree from an alignment

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/raxml:8.2.11_bcw_0.3.0'

    # Input of an alignment in FASTA format
    in_align_fasta = None

    # Parameter: Path + filename where the resultant tree should go
    tree_path = sl.Parameter()
    tree_stats_path = sl.Parameter()
    # DIRECTORY where the intermediate RAxML files should go (container fs space)
    raxml_working_dir = sl.Parameter(default=os.path.join(
        '/tmp',
        str(uuid.uuid4())
    ))

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
                    ' && cp $raxml_working_dir/RAxML_info.%s $out_tree_stats' % name +
                    ' && rm -r $raxml_working_dir',
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
            try:
                os.makedirs(os.path.dirname(self.out_accessions().path))
            except FileExistsError:
                pass
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

    chunk_size = sl.Parameter(default=40)

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
            gc.collect()


class NT_Repo_Prokka(sl.ContainerTask):
    container = 'golob/metaspades:v3.11.1--8A__bcw__0.3.0'
    in_repo = None
    workdir = sl.Parameter()
    num_concurrent = sl.Parameter(1)

    def out_repo(self):
        return NCBI_Repo_Peptides_TargetInfo(
            self,
            self.in_repo().path,
            )

    def annotate_version(self, version):
        try:
            seq_T = sl.ContainerTargetInfo(
                self,
                os.path.join(
                    self.workdir,
                    version,
                    'seq.fasta'
                ),
            )
            faa_T = sl.ContainerTargetInfo(
                self,
                os.path.join(
                    self.workdir,
                    version,
                    'peptides.faa'
                ),
            )
            ffn_T = sl.ContainerTargetInfo(
                self,
                os.path.join(
                    self.workdir,
                    version,
                    'rrna.ffn'
                ),
            )
            with seq_T.open('w') as seq_h:
                seq_h.write(
                        self.in_repo().target.repo.get_full_sequence(version)
                    )

            input_targets = {
                'seq': seq_T,
            }
            output_targets = {
                'faa': faa_T,
                'ffn': ffn_T,
            }

            self.ex(
                command=(
                    'prokka '
                    '--noanno '
                    '--outdir /share/$version '
                    '--prefix prokka '
                    '--cpus $vcpu '
                    '--force '
                    '$seq && '
                    'mv /share/$version/prokka.faa $faa && '
                    'mv /share/$version/prokka.ffn $ffn && '
                    'rm -r /share/$version '
                ),
                input_targets=input_targets,
                output_targets=output_targets,
                extra_params={
                    'version': version,
                    'vcpu': self.containerinfo.vcpu
                }
            )
            with faa_T.open() as faa_h:
                peptides = list({str(sr.seq) for sr in SeqIO.parse(faa_h, 'fasta')})
            with ffn_T.open() as ffn_h:
                rrna_16s = list({
                    str(sr.seq) for sr in SeqIO.parse(ffn_h, 'fasta')
                    if '16S' in sr.description
                    })
            log.info("Adding {} peptides to {}".format(
                len(peptides),
                version
            ))
            self.out_repo().target.repo.add_peptides_to_version(
                version,
                peptides
            )
            log.info("Adding {} 16S rRNA to {}".format(
                len(rrna_16s),
                version
            ))
            self.out_repo().target.repo.add_rRNA16s_to_version(
                version,
                rrna_16s
            )
            faa_T.target.remove()
            ffn_T.target.remove()
            seq_T.target.remove()
        except Exception as e:
            log.error(e)


    def run(self):
        versions_to_annotate = list(self.in_repo().target.repo.versions_needing_peptides())
        log.info("Prokka annotating {} genomes".format(len(versions_to_annotate)))
        annotate_pool = ThreadPool(
            int(self.num_concurrent))
        annotate_pool.map(
            func=self.annotate_version,
            iterable=versions_to_annotate,
        )


class NT_Repo_Output_FastaSeqInfo(sl.Task):
    in_repo = None

    fn_fasta_gz = sl.Parameter()
    fn_seq_info = sl.Parameter()
    min_len = sl.Parameter(default=1200)
    max_ambiguous_pct = sl.Parameter(default=1)

    def out_seqs(self):
        return sl.ContainerTargetInfo(
            self,
            self.fn_fasta_gz,
            format=luigi.format.Gzip
        )

    def out_seq_info(self):
        return sl.ContainerTargetInfo(
            self,
            self.fn_seq_info,
        )

    def run(self):
        repo = self.in_repo().target.repo

        with self.out_seqs().open('w') as out_h:
            with self.out_seq_info().open('w') as seq_info_h:
                seq_info_w = csv.DictWriter(
                    seq_info_h,
                    fieldnames=[
                        'seqname',
                        'version',
                        'accession',
                        'name',
                        'description',
                        'gi',
                        'tax_id',
                        'date',
                        'source',
                        'keywords',
                        'organism',
                        'length',
                        'ambig_count',
                        'seq_start',
                        'seq_stop',
                        'is_type',
                        'is_genome',
                    ])
                seq_info_w.writeheader()
                for i, rRNA_16s in enumerate(repo.find_feature(
                    qual_val='16S ribosomal RNA',
                    extra_fields=[
                        'tax_id',
                        'source',
                        'description',
                        'download_date',
                        'refseq__accession'])):
                        if i % 1000 == 0:
                            log.info("Exporting 16S rRNA {:,}".format(i))
                        if len(rRNA_16s['seq']) < self.min_len:
                            continue
                        # Implicit else
                        ambig_count = len([b for b in rRNA_16s['seq'].upper()
                                           if b not in {'A', 'C', 'G', 'T'}])
                        if (100.0 * ambig_count / len(rRNA_16s['seq']) > self.max_ambiguous_pct):
                            continue
                        # Implicit else basic test passed
                        name = "{}__{}__{}".format(
                            rRNA_16s['accession_version'].split('.')[0],
                            rRNA_16s['start'],
                            rRNA_16s['stop'])
                        out_h.write(">{} accession_version='{}'\n{}\n".format(
                            name,
                            rRNA_16s['accession_version'],
                            rRNA_16s['seq']
                        ).encode())
                        seq_info_w.writerow({
                            'seqname': name,
                            'version': rRNA_16s['accession_version'],
                            'accession': rRNA_16s['accession_version'].split('.')[0],
                            'name': rRNA_16s['accession_version'].split('.')[0],
                            'description': rRNA_16s['description'],
                            'gi': rRNA_16s['refseq__accession'],
                            'tax_id': rRNA_16s['tax_id'],
                            'date': rRNA_16s['download_date'],
                            'source': rRNA_16s['source'],
                            'keywords': "",
                            'organism': rRNA_16s['source'],
                            'length': len(rRNA_16s['seq']),
                            'ambig_count': ambig_count,
                            'seq_start': rRNA_16s['start'],
                            'seq_stop': rRNA_16s['stop'],
                            'is_type': "",
                            'is_genome': True,
                        })


class BuildTaxtasticDB(sl.ContainerTask):

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'

    # Where to put the sqlite database
    tax_db_path = sl.Parameter()

    def out_tax_db(self):
        return sl.ContainerTargetInfo(self, self.tax_db_path)

    def run(self):
        # Get our host paths for inputs and outputs
        output_targets = {
            'tax_db': self.out_tax_db(),
        }
        download_dir = os.path.join(
            self.containerinfo.container_working_dir,
            str(uuid.uuid4())
        )

        self.ex(
            command='mkdir -p $download_dir && '
                    ' taxit ' +
                    ' new_database' +
                    ' $tax_db'
                    ' -p $download_dir '
                    ' && rm -r $download_dir'
                    ,
            output_targets=output_targets,
            extra_params={
                'download_dir': download_dir,
            })


class FilterSeqinfoToFASTA(sl.Task):
    # Filter a sequence info csv to fit only sequences in a fasta file
    in_fasta = None
    in_seq_info = None

    filtered_seq_info_fn = sl.Parameter()

    def out_seq_info(self):
        return sl.ContainerTargetInfo(self, self.filtered_seq_info_fn)

    def run(self):
        # Get the sequence IDs in the fasta file
        seq_ids = {sr.id for sr in SeqIO.parse(self.in_fasta().open(), 'fasta')}

        with self.in_seq_info().open('r') as seq_info_h:
            seq_info_reader = csv.DictReader(seq_info_h)
            with self.out_seq_info().open('w') as filtered_seq_info_h:
                filtered_seq_info_w = csv.DictWriter(
                    filtered_seq_info_h,
                    fieldnames=seq_info_reader.fieldnames)
                filtered_seq_info_w.writeheader()
                for row in seq_info_reader:
                    if row['seqname'] in seq_ids:
                        filtered_seq_info_w.writerow(row)


class TaxTableForSeqInfo(sl.ContainerTask):
    # A Task that uses vsearch to find matches for experimental sequences in a repo of sequences

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'

    taxtable_path = sl.Parameter()

    in_tax_db = None
    in_seq_info = None

    def out_taxtable(self):
        return sl.ContainerTargetInfo(self, self.taxtable_path)

    def run(self):
        # Get our host paths for inputs and outputs
        input_targets = {
            'tax_db': self.in_tax_db(),
            'seq_info': self.in_seq_info()
        }
        output_targets = {
            'taxtable': self.out_taxtable(),
        }

        self.ex(
            command='taxit ' +
                    ' taxtable' +
                    ' $tax_db' +
                    ' --seq-info $seq_info' +
                    ' --outfile $taxtable',
            input_targets=input_targets,
            output_targets=output_targets,
            )


class ObtainCM(sl.ContainerTask):
    # Grab the cm alignment file from the cmalign container
    container = 'golob/infernal:1.1.2_bcw_0.2.0'

    cm_destination = sl.Parameter()

    def out_cm(self):
        return sl.ContainerTargetInfo(self, self.cm_destination)

    def run(self):
        output_targets = {
            'cm_dest': self.out_cm()
        }

        self.ex(
            command='cp /cmalign/data/SSU_rRNA_bacteria.cm $cm_dest',
            output_targets=output_targets
        )


class CombineRefpkg(sl.ContainerTask):
    # A task to take all the components of a refpkg and complete with a contents JSON
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0D'

    # Parameters to tell us where the refpkg should be placed,
    # as in: refpkg_path/refpkg_name.refpkg_version.refpkg
    refpkg_path = sl.Parameter()
    # what should we call this refpkg
    refpkg_name = sl.Parameter()
    # version, defaults to a timestamp of YYYYmmdd_HHMMSS
    refpkg_version = sl.Parameter(default=datetime.now().strftime('%Y%m%d_%H%M%S'))

    working_dir = sl.Parameter(default=os.path.join(
        '/tmp',
        str(uuid.uuid4())
    ))


    # dependencies
    in_aln_fasta = None
    in_aln_sto = None
    in_tree = None
    in_tree_stats = None
    in_taxtable = None
    in_seq_info = None
    in_cm = None

    def out_refpkg_tgz(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.refpkg_path,
                "{}.tgz".format(self.refpkg_name_ver())
            ),
            format=luigi.format.Gzip,
        )

    def refpkg_name_ver(self):
        return "{}.{}.refpkg".format(self.refpkg_name, self.refpkg_version)

    def run(self):
        input_targets = {
            'aln_fasta': self.in_aln_fasta(),
            'aln_sto': self.in_aln_sto(),
            'tree': self.in_tree(),
            'tree_stats': self.in_tree_stats(),
            'taxtable': self.in_taxtable(),
            'seq_info': self.in_seq_info(),
            'cm': self.in_cm(),
        }
        output_targets = {
            'refpkg_tgz': self.out_refpkg_tgz()
        }
        # how to handle directory target?

        self.ex(
            command='mkdir -p $working_dir && cd $working_dir && taxit' +
            ' create --locus 16S' +
            ' --package-name $refpkg_name --clobber' +
            ' --aln-fasta $aln_fasta ' +
            ' --aln-sto $aln_sto ' +
            ' --tree-file $tree ' +
            ' --tree-stats $tree_stats ' +
            ' --taxonomy $taxtable ' +
            ' --seq-info $seq_info ' +
            ' --profile $cm ' +
            ' && tar czvf $refpkg_tgz $refpkg_name/* '+
            ' && rm -r $working_dir ',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'refpkg_name': self.refpkg_name_ver(),
                'working_dir': self.working_dir,
            }
        )


class ExtractRefpkgAlignment(sl.Task):
    # When given a refpkg in tar.gz format, unpack to the relevant outputs

    in_refpkg_tgz = None

    # Where to extract our alignments
    aln_fasta_fn = sl.Parameter()
    aln_sto_fn = sl.Parameter()

    def out_aln_fasta(self):
        return sl.ContainerTargetInfo(
            self,
            self.aln_fasta_fn
        )

    def out_aln_sto(self):
        return sl.ContainerTargetInfo(
            self,
            self.aln_sto_fn
        )

    def run(self):
        # Open the tarfile
        with self.in_refpkg_tgz().target.open('rb') as refpkg_tgz_h:
            tar_h = tarfile.open(
                mode='r:*',
                fileobj=io.BytesIO(refpkg_tgz_h.read())
            )
            # Get the contents keyed by the filenames, excluding dirs
            tar_contents_dict = {os.path.basename(f.name): f for f in tar_h.getmembers()}
            try:
                contents = json.load(tar_h.extractfile(tar_contents_dict['CONTENTS.json']))
            except TypeError:
                contents = json.loads(
                    tar_h.extractfile(
                            tar_contents_dict['CONTENTS.json']
                        ).read().decode('utf-8')
                    )
            aln_fasta_intgz = contents['files'].get('aln_fasta')
            aln_sto_intgz = contents['files'].get('aln_sto')

            if aln_fasta_intgz and aln_sto_intgz:
                # Both version of the alignment are in the refpkg
                with self.out_aln_fasta().open('w') as out_aln_fasta_h:
                    out_aln_fasta_h.write(
                        tar_h.extractfile(
                            tar_contents_dict[aln_fasta_intgz]
                        ).read().decode('utf-8')
                    )
                with self.out_aln_sto().open('w') as out_aln_sto_h:
                    out_aln_sto_h.write(
                        tar_h.extractfile(
                            tar_contents_dict[aln_sto_intgz]
                        ).read().decode('utf-8')
                    )
            elif aln_fasta_intgz:
                # Only fasta exists
                with self.out_aln_fasta().open('w') as out_aln_fasta_h:
                    out_aln_fasta_h.write(
                        tar_h.extractfile(
                            tar_contents_dict[aln_fasta_intgz]
                        ).read().decode('utf-8')
                    )
                # And convert to sto format
                with self.out_aln_sto().open('w') as out_aln_sto_h:
                    AlignIO.write(
                        AlignIO.read(
                            tar_h.extractfile(tar_contents_dict[aln_fasta_intgz]),
                            'fasta'),
                        out_aln_sto_h,
                        'stockholm'
                    )
            elif aln_sto_intgz:
                # Only STO exists
                with self.out_aln_sto().open('w') as out_aln_sto_h:
                    out_aln_sto_h.write(
                        tar_h.extractfile(
                            tar_contents_dict[aln_sto_intgz]
                        ).read().decode('utf-8')
                    )
                with self.out_aln_fasta().open('w') as out_aln_fasta_h:
                    AlignIO.write(
                        AlignIO.read(
                                     tar_h.extractfile(tar_contents_dict[aln_sto_intgz]),
                                     'stockholm'),
                        out_aln_fasta_h,
                        'fasta'
                    )
            else:
                # NO alignment present
                raise Exception("Refset at {} does not contain an alignment".format(
                    self.in_refpkg_tgz
                ))


class CombineAlignmentsSTO(sl.ContainerTask):
    # Grab the cm alignment file from the cmalign container
    container = 'golob/infernal:1.1.2_bcw_0.2.0'

    in_aln_sto_1 = None
    in_aln_sto_2 = None

    combined_aln_sto_fn = sl.Parameter()

    def out_aln_sto(self):
        return sl.ContainerTargetInfo(self, self.combined_aln_sto_fn)

    def run(self):
        input_targets = {
            'in_aln_1': self.in_aln_sto_1(),
            'in_aln_2': self.in_aln_sto_2()
        }
        output_targets = {
            'out_aln': self.out_aln_sto()
        }

        self.ex(
            command='esl-alimerge --dna' +
                    ' -o $out_aln' +
                    ' $in_aln_1 $in_aln_2',
            input_targets=input_targets,
            output_targets=output_targets
        )


class PPLACER_PlaceAlignment(sl.ContainerTask):
    # Place aligned refpkg and sv onto a refpkg tree
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    #  Where to put the resultant jplace output
    jplace_fn = sl.Parameter()
    # pplacer parameters
    pplacer_prior_lower = sl.Parameter(default=0.01)

    #  Dependencies
    in_refpkg_tgz = None
    in_merged_aln_sto = None

    def out_jplace(self):
        return sl.ContainerTargetInfo(self, self.jplace_fn)

    def run(self):
        input_targets = {
            'refpkg_tgz': self.in_refpkg_tgz(),
            'merged_aln': self.in_merged_aln_sto(),
        }
        output_targets = {
            'jplace': self.out_jplace()
        }

        refpkg_rel_path = self.in_refpkg_tgz().get_refpkg_rel_path()

        unique_prefix = str(uuid.uuid4())
        run_tmpdir = "/tmp/{}".format(unique_prefix)
        self.ex(
            command='mkdir -p $run_tmpdir/refpkg/ && tar xzvf $refpkg_tgz -C $run_tmpdir/refpkg/ ' +
                    ' && pplacer -p --inform-prior --prior-lower $prior_lower --map-identity' +
                    ' -j $nproc' +
                    ' -c $run_tmpdir/refpkg/$refpkg_rel_path' +
                    ' $merged_aln' +
                    ' -o $jplace' +
                    ' && rm -r $run_tmpdir',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'prior_lower': self.pplacer_prior_lower,
                'nproc': self.containerinfo.vcpu,
                'refpkg_rel_path': refpkg_rel_path,
                'run_tmpdir': run_tmpdir,
            }
        )


class Jplace_Reduplicate(sl.ContainerTask):
    # Reduplicate a jplace
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    # Parameters
    jplace_fn = sl.Parameter()

    # Dependencies
    in_jplace = None
    in_weights = None

    def out_jplace(self):
        # Gzip format by default
        return sl.ContainerTargetInfo(
                                self,
                                self.jplace_fn,
                                format=luigi.format.Gzip,
        )

    def run(self):
        input_targets = {
            'in_jplace': self.in_jplace(),
            'weights': self.in_weights()
        }
        output_targets = {
            'out_jplace': self.out_jplace()
        }

        self.ex(
            command='guppy redup -m' +
                    ' -o /dev/stdout' +
                    ' -d $weights' +
                    ' $in_jplace' +
                    ' | gzip > $out_jplace',
            input_targets=input_targets,
            output_targets=output_targets
        )


class Jplace_PCA(sl.ContainerTask):
    # Calculate EPCA for a JPLACE
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_jplace = None
    in_refpkg_tgz = None
    in_seq_map = None

    # Parameter for prefix / path
    # guppy will take /path/to/prefix
    # and create /path/to/prefix.proj
    #           /path/to/prefix.trans
    #           /path/to/prefix.xml
    path = sl.Parameter()
    pca = sl.Parameter()
    prefix = sl.Parameter()

    def out_proj(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.proj".format(self.prefix)
            )
        )

    def out_trans(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.trans".format(self.prefix)
            )
        )

    def out_xml(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.xml".format(self.prefix)
            )
        )

    def run(self):
        if not (self.pca.lower() == 'lpca' or self.pca.lower() == 'epca'):
            raise Exception("{} is not a valid PCA type (EPCA or LPCA)".format(self.pca))
        pca = self.pca.lower()

        refpkg_rel_path = self.in_refpkg_tgz().get_refpkg_rel_path()
        input_targets = {
            'jplace': self.in_jplace(),
            'seq_map': self.in_seq_map(),
            'refpkg_tgz': self.in_refpkg_tgz(),
        }
        output_targets = {
            'proj': self.out_proj(),
            'trans': self.out_trans(),
            'xml': self.out_xml(),
        }

        unique_prefix = str(uuid.uuid4())
        run_tmpdir = "/tmp/{}".format(unique_prefix)
        self.ex(
            command='mkdir -p $run_tmpdir/refpkg/ && cd $run_tmpdir/refpkg/ && tar xzvf $refpkg_tgz -C $run_tmpdir/refpkg/ ' +
                    ' && mkdir -p $run_tmpdir/results ' +
                    ' && guppy $pca' +
                    ' $jplace:$seq_map ' +
                    ' -c $run_tmpdir/refpkg/$refpkg_rel_path ' + 
                    ' --out-dir $run_tmpdir/results --prefix $prefix' +
                    ' && cp $run_tmpdir/results/$prefix.xml $xml' + 
                    ' && cp $run_tmpdir/results/$prefix.proj $proj' + 
                    ' && cp $run_tmpdir/results/$prefix.trans $trans' +
                    ' && rm -r $run_tmpdir',

            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'refpkg_rel_path': refpkg_rel_path,
                'prefix': self.prefix,
                'pca': pca,
                'run_tmpdir': run_tmpdir
            }
        )


class Jplace_ADCL(sl.ContainerTask):
    # Calculate ADCL for a JPLACE
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_jplace = None

    adcl_fn = sl.Parameter()

    def out_adcl_gz(self):
        return sl.ContainerTargetInfo(
                                self,
                                self.adcl_fn,
                                format=luigi.format.Gzip,
        )

    def run(self):
        input_targets = {
            'in_jplace': self.in_jplace()
        }

        output_targets = {
            'adcl_gz': self.out_adcl_gz()
        }

        self.ex(
            command='(echo name,adcl,weight && '
                    ' guppy adcl --no-collapse $in_jplace -o /dev/stdout) | '
                    ' gzip > $adcl_gz',
            input_targets=input_targets,
            output_targets=output_targets
        )


class Jplace_EDPL(sl.ContainerTask):
    # Calculate ADCL for a JPLACE
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_jplace = None

    edpl_fn = sl.Parameter()

    def out_edpl_gz(self):
        return sl.ContainerTargetInfo(
                                self,
                                self.edpl_fn,
                                format=luigi.format.Gzip,
        )

    def run(self):
        input_targets = {
            'in_jplace': self.in_jplace()
        }

        output_targets = {
            'edpl_gz': self.out_edpl_gz()
        }

        self.ex(
            command='(echo name,edpl && '
                    ' guppy edpl --csv $in_jplace -o /dev/stdout) | '
                    ' gzip > $edpl_gz',
            input_targets=input_targets,
            output_targets=output_targets
        )


class Jplace_KR_Distance(sl.ContainerTask):
    # Calculate KR phylogenetic distance for a placement
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_jplace = None
    in_refpkg_tgz = None
    in_seq_map = None

    kr_fn = sl.Parameter()

    def out_kr_distance(self):
        return sl.ContainerTargetInfo(
            self,
            self.kr_fn
        )

    def run(self):
        refpkg_rel_path = self.in_refpkg_tgz().get_refpkg_rel_path()
        input_targets = {
            'jplace': self.in_jplace(),
            'seq_map': self.in_seq_map(),
            'refpkg_tgz': self.in_refpkg_tgz(),
        }
        output_targets = {
            'kr_distance': self.out_kr_distance(),
        }

        unique_prefix = str(uuid.uuid4())
        run_tmpdir = "/tmp/{}".format(unique_prefix)
        self.ex(
            command='mkdir -p $run_tmpdir/refpkg/ && cd $run_tmpdir/refpkg/ && tar xzvf $refpkg_tgz -C $run_tmpdir/refpkg/ ' +
                    ' && guppy kr' +
                    ' --list-out' +
                    ' -c $run_tmpdir/refpkg/$refpkg_rel_path ' +
                    ' $jplace:$seq_map ' +
                    ' -o $kr_distance' +
                    ' && rm -r $run_tmpdir',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'refpkg_rel_path': refpkg_rel_path,
                'run_tmpdir': run_tmpdir,
            }
        )


class Jplace_Alpha_Diversity(sl.ContainerTask):
    # Calculate alpha diversity for a placement
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_jplace = None
    in_seq_map = None

    alpha_diversity_fn = sl.Parameter()

    def out_alpha_diversity(self):
        return sl.ContainerTargetInfo(
            self,
            self.alpha_diversity_fn
        )

    def run(self):
        input_targets = {
            'jplace': self.in_jplace(),
            'seq_map': self.in_seq_map(),
        }
        output_targets = {
            'alpha_diversity': self.out_alpha_diversity(),
        }

        self.ex(
            command='guppy fpd' +
                    ' --csv' +
                    ' --include-pendant' +
                    ' --chao-d 0,1,1.00001,2,3,4,5 '
                    ' $jplace:$seq_map ' +
                    ' -o $alpha_diversity',
            input_targets=input_targets,
            output_targets=output_targets,
        )


class PlacementDB_Prep(sl.ContainerTask):
    # Prep a placement db
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_refpkg_tgz = None

    placement_db_fn = sl.Parameter()

    def out_placement_db(self):
        tax_ids = self.in_refpkg_tgz().get_tax_ids()
        ranks = self.in_refpkg_tgz().get_tax_ranks()
        return PlacementDB_Prepped_ContainerTargetInfo(
            self,
            self.placement_db_fn,
            expected_tax_ids=tax_ids,
            expected_ranks=ranks,
        )

    def run(self):
        refpkg_rel_path = self.in_refpkg_tgz().get_refpkg_rel_path()
        input_targets = {
            'refpkg_tgz': self.in_refpkg_tgz(),
        }

        output_targets = {
            'placement_db': self.out_placement_db(),
        }

        unique_prefix = str(uuid.uuid4())
        run_tmpdir = "/tmp/{}".format(unique_prefix)
        self.ex(
            command='mkdir -p $run_tmpdir/refpkg/ && tar xzvf $refpkg_tgz -C $run_tmpdir/refpkg/ ' +
                    ' && rppr prep_db' +
                    ' -c $run_tmpdir/refpkg/$refpkg_rel_path ' +
                    ' --sqlite $placement_db' +
                    '&& rm -r $run_tmpdir ',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'refpkg_rel_path': refpkg_rel_path,
                'run_tmpdir': run_tmpdir
            }
        )


class PlacementDB_Classify_SV(sl.ContainerTask):
    # Prep a placement db
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    # Dependencies
    in_refpkg_tgz = None
    in_sv_refpkg_aln_sto = None
    in_placement_db = None
    in_jplace = None

    # Parameters
    classifer = sl.Parameter(default='hybrid2')
    seed = sl.Parameter(default='1')
    likelihood_cutoff = sl.Parameter(default='0.9')
    bayes_cutoff = sl.Parameter(default='1.0')
    multiclass_min = sl.Parameter(default='0.2')
    bootstrap_cutoff = sl.Parameter(default='0.8')
    bootstrap_extension_cutoff = sl.Parameter(default='0.4')
    nbc_word_length = sl.Parameter(default='8')
    nbc_target_rank = sl.Parameter(default='genus')
    nbc_boot = sl.Parameter(default='100')
    use_posterior_prob = sl.Parameter(default=True)
    no_pre_mask = sl.Parameter(default=False)
    no_random_tie_break = sl.Parameter(default=False)

    def out_placement_db(self):
        return PlacementDB_Classified_ContainerTargetInfo(
            self,
            self.in_placement_db().path,
            guppy_command=self.make_guppy_command(
                    self.guppy_parameter_to_dict()
                )
        )

    def guppy_parameter_to_dict(self):
        return {
            'classifier': self.classifer,
            'nproc': self.containerinfo.vcpu,
            'seed': self.seed,
            'likelihood_cutoff': self.likelihood_cutoff,
            'bayes_cutoff': self.bayes_cutoff,
            'multiclass_min': self.multiclass_min,
            'bootstrap_cutoff': self.bootstrap_cutoff,
            'bootstrap_extension_cutoff': self.bootstrap_extension_cutoff,
            'nbc_word_length': self.nbc_word_length,
            'nbc_target_rank': self.nbc_target_rank,
            'nbc_boot': self.nbc_boot,
            }

    def make_guppy_command(self, replace_dict=None):
        s = (
            ' guppy classify'
            ' --classifier $classifier '
            ' -j $nproc '
            ' -c $run_tmpdir/refpkg/$refpkg_rel_path '
            ' --nbc-sequences $sv_refpkg_aln_sto '
            ' --sqlite $placement_db'
            ' --seed $seed'
            ' --cutoff $likelihood_cutoff'
            ' --bayes-cutoff $bayes_cutoff'
            ' --multiclass-min $multiclass_min'
            ' --bootstrap-cutoff $bootstrap_cutoff'
            ' --bootstrap-extension-cutoff $bootstrap_extension_cutoff'
            ' --word-length $nbc_word_length'
            ' --nbc-rank $nbc_target_rank'
            ' --n-boot $nbc_boot'
            ' $jplace'
        )
        if self.use_posterior_prob:
            s += ' --pp'
        if self.no_pre_mask:
            s += ' --no-pre-mask'
        if self.no_random_tie_break:
            s += ' --no-random-tie-break'

        if replace_dict:
            return Template(s).safe_substitute(replace_dict)
        else:
            return s

    def run(self):
        unique_prefix = str(uuid.uuid4())
        run_tmpdir = "/tmp/{}".format(unique_prefix)
        refpkg_rel_path = self.in_refpkg_tgz().get_refpkg_rel_path()
        input_targets = {
            'refpkg_tgz': self.in_refpkg_tgz(),
            'sv_refpkg_aln_sto': self.in_sv_refpkg_aln_sto(),
            'jplace': self.in_jplace(),
        }

        output_targets = {
            'placement_db': self.out_placement_db(),
        }

        extra_params = {
                'refpkg_rel_path': refpkg_rel_path,
                'run_tmpdir': run_tmpdir,
            }
        extra_params.update(self.guppy_parameter_to_dict())

        self.ex(
            command='mkdir -p $run_tmpdir/refpkg/ && tar xzvf $refpkg_tgz -C $run_tmpdir/refpkg/ &&' +
            self.make_guppy_command()
            + ' && rm -r $run_tmpdir ',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params=extra_params,
        )


class PlacementDB_MCC(sl.ContainerTask):
    # Prep a placement db
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_placement_db = None
    in_weights = None

    def out_placement_db(self):
        return PlacementDB_MCC_ContainerTargetInfo(
            self,
            self.in_placement_db().path,
            expected_seq_ids=self.all_seq_ids()
        )

    def all_seq_ids(self):
        with self.in_weights().open() as seq_weights_h:
            return {
                r[1] for r in csv.reader(seq_weights_h)
            }

    def run(self):  
        input_targets = {
            'weights': self.in_weights()
        }

        output_targets = {
            'placement_db': self.out_placement_db(),
        }

        self.ex(
            command=(
                'multiclass_concat.py '
                ' -k'
                ' --dedup-info $weights'
                ' $placement_db'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
        )


class PlacementDB_AddSI(sl.ContainerTask):
    # Prep a placement db
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_placement_db = None
    in_seq_map = None

    def out_placement_db(self):
        return PlacementDB_SI_ContainerTargetInfo(
            self,
            self.in_placement_db().path,
            expected_seq_ids=self.all_seq_ids(),
        )

    def all_seq_ids(self):
        with self.in_seq_map().open() as seq_map_h:
            return {
                r[0] for r in csv.reader(seq_map_h)
            }

    def run(self):
        input_targets = {
            'seq_map': self.in_seq_map()
        }

        output_targets = {
            'placement_db': self.out_placement_db(),
        }

        self.ex(
            command=(
                '(echo "name,specimen"; cat $seq_map) | '
                ' csvsql --table seq_info '
                ' --insert --snifflimit 1000'
                ' --db sqlite:///$placement_db'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
        )


class GenerateTables(sl.ContainerTask):
    # by specimen tabular output
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0c'

    in_placement_db = None
    in_seq_map = None
    in_labels = None

    tables_path = sl.Parameter()
    rank = sl.Parameter()

    def out_by_specimen(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.tables_path,
                'by_specimen.{}.csv'.format(self.rank)
            )
        )

    def out_by_taxon(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.tables_path,
                'by_taxon.{}.csv'.format(self.rank)
            )
        )

    def out_tallies_wide(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.tables_path,
                'tallies_wide.{}.csv'.format(self.rank)
            )
        )

    def run(self):
        unique_prefix = str(uuid.uuid4())
        run_tmpdir = "/tmp/{}".format(unique_prefix)

        input_targets = {
            'placement_db': self.in_placement_db(),
            'seq_map': self.in_seq_map(),
        }

        if self.in_labels:
            input_targets['labels'] = self.in_labels()

        output_targets = {
            'by_specimen': self.out_by_specimen(),
            'by_taxon': self.out_by_taxon(),
            'tallies_wide': self.out_tallies_wide(),
        }

        command = (
                'mkdir -p $run_tmpdir/tables &&'
                ' python2 /usr/bin/classif_table.py'
                ' $placement_db'
                ' $run_tmpdir/tables/by_taxon.csv'
                ' --rank $rank'
                ' --specimen-map $seq_map'
                ' --by-specimen $run_tmpdir/tables/by_specimen.csv'
                ' --tallies-wide $run_tmpdir/tables/tallies_wide.csv'
                ' && cp -r $run_tmpdir/tables/by_taxon.csv $by_taxon'
                ' && cp -r $run_tmpdir/tables/by_specimen.csv $by_specimen'
                ' && cp -r $run_tmpdir/tables/tallies_wide.csv $tallies_wide'
        )

        if self.in_labels:
            command += ' --metadata-map $labels'
        command += ' && rm -r $run_tmpdir/'
        self.ex(
            command=command,
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'rank': self.rank,
                'run_tmpdir': run_tmpdir,
            }
        )


class LoadManifest(LoadFile):

    def manifest(self):
        # Load manifest as csv and return as a list of dicts
        with self.out_file().open() as manifest_h:
            return [r for r in csv.DictReader(manifest_h)]

    def batched_specimens(self):
        manifest = self.manifest()
        batches = {r.get('batch') for r in manifest}
        log.info("{} batches".format(
            len(batches))
            )

        for batch in batches:
            yield((batch, {
                r['specimen'] for r in manifest
                if r.get('batch') == batch
            }))

    def get_columns(self):
        with self.out_file().open() as manifest_h:
            return set(csv.DictReader(manifest_h).fieldnames)

    def is_paired(self):
        return 'read__2' in self.get_columns()

    def has_index(self):
        columns = self.get_columns()
        if 'read__2' in columns and 'index__2' in columns and 'index__1' in columns:
            return True
        elif 'index__1' in columns:
            return True
        else:
            return False

    def is_valid(self):
        columns = self.get_columns()
        if 'read__1' in columns and 'specmen' in columns:
            return True
        else:
            return False

    def has_batches(self):
        return 'batch' in self.get_columns()

    def get_specimens(self):
        with self.out_file().open() as manifest_h:
            return {r.get('specimen') for r in csv.DictReader(manifest_h)}


class LoadSpecimenReads(sl.ExternalTask):
    # Given a manifest and a specimen ID to target, load the read(s)
    in_manifest = None
    specimen = sl.Parameter()

    def specimen_manifest(self):
        with self.in_manifest().open() as manifest_h:
            for r in csv.DictReader(manifest_h):
                if r.get('specimen') == self.specimen:
                    return r
        # Implicit else we didn't find the specimen...
        raise Exception("Could not find specimen {} in the manifest".format(
            self.specimen
        ))

    def out_reads(self):
        reads_dict = {}
        specimen_manifest = self.specimen_manifest()
        if not 'read__1' in specimen_manifest:
            raise Exception("No forward / primary read for specimen {}".format(
                self.specimen
            ))
        reads_dict['R1'] = sl.ContainerTargetInfo(
            self,
            specimen_manifest['read__1'],
            format=luigi.format.Nop
        )
        if 'read__2' in specimen_manifest:
            reads_dict['R2'] = sl.ContainerTargetInfo(
                self,
                specimen_manifest['read__2'],
                format=luigi.format.Nop
            )
        if 'index__1' in specimen_manifest and specimen_manifest['index__1'] != "":
            reads_dict['I1'] = sl.ContainerTargetInfo(
                self,
                specimen_manifest['index__1'],
                format=luigi.format.Nop
            )
        if 'index__2' in specimen_manifest and specimen_manifest['index__2'] != "":
            reads_dict['I2'] = sl.ContainerTargetInfo(
                self,
                specimen_manifest['index__2'],
                format=luigi.format.Nop
            )
        return reads_dict


class BCCSpecimenReads(sl.ContainerTask):
    # Use barcodecop to verify reads were properly demultiplexed
    container = 'golob/barcodecop:0.4.1__bcw_0.3.0'

    # For dependencies
    in_reads = None

    specimen = sl.Parameter()
    path = sl.Parameter()

    def out_reads(self):
        reads_dict = {}
        reads_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.bcc.fq.gz".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        reads_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.bcc.fq.gz".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        return reads_dict
    """
    def out_stats(self):
        stats_dict = {}
        stats_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.bcc.stats.csv".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        stats_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.bcc.stats.csv".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        return stats_dict
    """

    def run(self):
        input_targets={
            'read_1': self.in_reads().get('R1'),
            'read_2': self.in_reads().get('R2'),
            'index_1': self.in_reads().get('I1'),
            'index_2': self.in_reads().get('I2'),
        }
        output_targets={
            'bcc_read_1': self.out_reads()['R1'],
            'bcc_read_2': self.out_reads()['R2'],
            #'bcc_stat_1': self.out_stats()['R1'],
            #'bcc_stat_2': self.out_stats()['R2']
        }
        self.ex(
            command=(
                'barcodecop'
                ' $index_1 $index_2'
                ' --match-filter'
                ' -f $read_1'
                ' -o $bcc_read_1'
            #    ' -C $bcc_stat_1'
                ' && barcodecop'
                ' $index_1 $index_2'
                ' --match-filter'
                ' -f $read_2'
                ' -o $bcc_read_2'
            #    ' -C $bcc_stat_2'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
        )

class DADA2_FilterAndTrim(sl.ContainerTask):
    # DADA2 filter and trim for a specimen (F and R)
    container = 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'

    # Dependencies
    in_reads = None

    # Parameters
    specimen = sl.Parameter()
    path = sl.Parameter()
    # DADA2 parameters
    maxN = sl.Parameter(default=0)
    maxEE = sl.Parameter(default='Inf')
    f_trunc = sl.Parameter(default=250)
    r_trunc = sl.Parameter(default=250)
    trim_left = sl.Parameter(default=15)
    trunc_Q = sl.Parameter(default=2)

    def out_reads(self):
        reads_dict = {}
        reads_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.ft.fq.gz".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        reads_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.ft.fq.gz".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        return reads_dict
    
    def run(self):
        input_targets = {
            'read_1': self.in_reads()['R1'],
            'read_2': self.in_reads()['R2']
        }
        output_targets = {
            'read_1_ft': self.out_reads()['R1'],
            'read_2_ft': self.out_reads()['R2']
        }

        self.ex(
            command=(
                'Rscript -e "'
                "library('dada2'); "
                'filterAndTrim('
                "'$read_1', '$read_1_ft', "
                "'$read_2', '$read_2_ft', "
                'trimLeft = $trim_left, '
                'maxN = $maxN, '
                'maxEE = c($maxEE, $maxEE), '
                'truncLen = c($f_trunc, $r_trunc), '
                'truncQ = $trunc_Q, '
                'compress = TRUE, '
                'verbose = TRUE, '
                'multithread = $vcpu)'
                '"'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'trim_left': self.trim_left,
                'maxN': self.maxN,
                'maxEE': self.maxEE,
                'f_trunc': self.f_trunc,
                'r_trunc': self.r_trunc,
                'trunc_Q': self.trunc_Q,
                'vcpu': self.containerinfo.vcpu
            }
        )


class DADA2_Dereplicate(sl.ContainerTask):
    # DADA2 filter and trim for a specimen (F and R)
    container = 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'

    # Dependencies
    in_reads = None

    # Parameters
    specimen = sl.Parameter()
    path = sl.Parameter()
    # DADA2 parameters
    maxN = sl.Parameter(default=1)
    maxEE = sl.Parameter(default='Inf')
    f_trunc = sl.Parameter(default=280)
    r_trunc = sl.Parameter(default=250)
    trim_left = sl.Parameter(default=15)
    trunc_Q = sl.Parameter(default=2)

    def out_rds(self):
        rds_dict = {}
        rds_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.derep.rds".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        rds_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.derep.rds".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        return rds_dict
    
    def run(self):
        input_targets = {
            'read_1': self.in_reads()['R1'],
            'read_2': self.in_reads()['R2']
        }
        output_targets = {
            'derep_1': self.out_rds()['R1'],
            'derep_2': self.out_rds()['R2']
        }

        self.ex(
            command=(
                'Rscript -e "'
                "library('dada2'); "
                "derep1 <- derepFastq('$read_1', verbose = TRUE); "
                "saveRDS(derep1, '$derep_1'); "
                "derep2 <- derepFastq('$read_2', verbose = TRUE); "
                "saveRDS(derep2, '$derep_2'); "
                '"'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
        )


class DADA2_LearnError(sl.ContainerTask):
    # DADA2 filter and trim for a specimen (F and R)
    container = 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'

    # Dependencies
    # Should be a LIST of DICTS 
    # of filter-trimmed read targets
    in_reads = None

    # Parameters
    batch = sl.Parameter()
    path = sl.Parameter()
    tar_reads = sl.Parameter(default=False)
    MAX_CONSIST = sl.Parameter(default=10)
    randomize = sl.Parameter(default='FALSE')
    n_bases = sl.Parameter(default='1e8')

    def out_rds(self):
        rds_dict = {}
        rds_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.err.rds".format(self.batch)
            ),
            format=luigi.format.Nop
        )
        rds_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.err.rds".format(self.batch)
            ),
            format=luigi.format.Nop
        )
        return rds_dict

    def out_csv(self):
        rds_dict = {}
        rds_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.err.csv".format(self.batch)
            ),
        )
        rds_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.err.csv".format(self.batch)
            ),
        )
        return rds_dict
    
    def run(self):
        if self.tar_reads == 'true':
            # Make a tarfile for the forward and reverse entries
            run_id = str(uuid.uuid4())
            F_reads_tar = sl.ContainerTargetInfo(
                self,
                os.path.join(
                    self.path,
                    "F_reads.{}.tar".format(self.batch)
                ),
                format=luigi.format.Nop
            )
            with F_reads_tar.open('w') as tar_h:
                arch = tarfile.open(fileobj=tar_h, mode='w:')
                for r_i, read_t in enumerate(self.in_reads):
                    with read_t()['R1'].open('r') as file_h:
                        try:
                            f_ti = arch.gettarinfo(
                                fileobj=file_h,
                                arcname='r_{}.fa.gz'.format(r_i)
                                )
                            arch.addfile(
                                f_ti,
                                fileobj=file_h
                            )
                        except AttributeError:  # S3 doesn't do well with this
                            # So manually download to a named-temp file and then add that to the tar
                            with tempfile.NamedTemporaryFile() as f_ntf:
                                f_ntf.write(file_h.read())
                                f_ntf.seek(0)
                                f_ti = arch.gettarinfo(
                                    fileobj=f_ntf,
                                    arcname='r_{}.fa.gz'.format(r_i)
                                )
                                arch.addfile(
                                    f_ti,
                                    fileobj=f_ntf
                                )

            input_targets_F = {
                'reads_tar': F_reads_tar,
            }

            R_reads_tar = sl.ContainerTargetInfo(
                self,
                os.path.join(
                    self.path,
                    "R_reads.{}.tar".format(self.batch)
                ),
                format=luigi.format.Nop
            )
            with R_reads_tar.open('w') as tar_h:
                arch = tarfile.open(fileobj=tar_h, mode='w:')
                for r_i, read_t in enumerate(self.in_reads):
                    with read_t()['R2'].open('r') as file_h:
                        try:
                            f_ti = arch.gettarinfo(
                                fileobj=file_h,
                                arcname='r_{}.fa.gz'.format(r_i)
                                )
                            arch.addfile(
                                f_ti,
                                fileobj=file_h
                            )
                        except AttributeError:  # S3 doesn't do well with this
                            # So manually download to a named-temp file and then add that to the tar
                            with tempfile.NamedTemporaryFile() as f_ntf:
                                f_ntf.write(file_h.read())
                                f_ntf.seek(0)
                                f_ti = arch.gettarinfo(
                                    fileobj=f_ntf,
                                    arcname='r_{}.fa.gz'.format(r_i)
                                )
                                arch.addfile(
                                    f_ti,
                                    fileobj=f_ntf
                                )

            input_targets_R = {
                'reads_tar': R_reads_tar,
            }

            output_targets_F = {
                'err_1': self.out_rds()['R1'],
                'err_csv_1': self.out_csv()['R1']
            }

            output_targets_R = {
                'err_2': self.out_rds()['R2'],
                'err_csv_2': self.out_csv()['R2'],
            }
            command = 'mkdir -p /tmp/{}/{} && '.format(run_id, self.batch) + \
                    'tar xf $reads_tar -C /tmp/{}/{}/ && '.format(run_id, self.batch) + \
                    'Rscript -e "' + \
                    "library('dada2'); "
            # Forward
            command_F = command + 'errF <- learnErrors(c('
            command_F += ",".join(["'/tmp/{}/{}/r_{}.fa.gz'".format(run_id, self.batch, i) for i in range(len(self.in_reads))])
            command_F += (
                '), multithread=$vcpu'
                "); saveRDS(errF, '$err_1'); "
                "write.csv(errF, '$err_csv_1');"
                '"'
            )
            # Reverse
            command_R = command + 'errR <- learnErrors(c('
            command_R += ",".join(["'/tmp/{}/{}/r_{}.fa.gz'".format(run_id, self.batch, i) for i in range(len(self.in_reads))])
            command_R += (
                '), multithread=$vcpu'
                "); saveRDS(errR, '$err_2'); "
                "write.csv(errR, '$err_csv_2');"
            )
            command_R += '"'
            # Forward:
            self.ex(
                command=command_F,
                input_targets=input_targets_F,
                output_targets=output_targets_F,
                extra_params={'vcpu': self.containerinfo.vcpu}
            )
            # Reverse:
            self.ex(
                command=command_R,
                input_targets=input_targets_R,
                output_targets=output_targets_R,
                extra_params={'vcpu': self.containerinfo.vcpu}
            )
            F_reads_tar.target.__del__()
            R_reads_tar.target.__del__()
        else: # Do not tar reads
            output_targets_F = {
                'err_1': self.out_rds()['R1'],
                'err_csv_1': self.out_csv()['R1']
            }
            input_targets_F = {
                'read_F_{}'.format(i): t()['R1'] for i, t in enumerate(self.in_reads)
            }
            output_targets_R = {
                'err_2': self.out_rds()['R2'],
                'err_csv_2': self.out_csv()['R2'],
            }
            input_targets_R = {
                'read_R_{}'.format(i): t()['R2'] for i, t in enumerate(self.in_reads)
            }


            command = (
                        'Rscript -e "'
                        "library('dada2'); "
            )

            # Forward
            command_F = command + 'errF <- learnErrors(c('
            command_F += ",".join(["'$read_F_{}'".format(i) for i in range(len(self.in_reads))])
            command_F += (
                '), multithread=$vcpu'
                ', MAX_CONSIST=$MAX_CONSIST'
                ', randomize=$randomize'
                ', nbases=$nbases'
                ', verbose=FALSE'
                "); saveRDS(errF, '$err_1'); "
                "write.csv(errF, '$err_csv_1');"
                '"'
            )
            # Reverse
            command_R = command + 'errR <- learnErrors(c('
            command_R += ",".join(["'$read_R_{}'".format(i) for i in range(len(self.in_reads))])
            command_R += (
                '), multithread=$vcpu'
                ', MAX_CONSIST=$MAX_CONSIST'
                ', randomize=$randomize'
                ', nbases=$nbases'
                ', verbose=FALSE'
                "); saveRDS(errR, '$err_2'); "
                "write.csv(errR, '$err_csv_2');"
            )
            command_R += '"'
            # Forward:
            self.ex(
                command=command_F,
                input_targets=input_targets_F,
                output_targets=output_targets_F,
                extra_params={
                    'vcpu': self.containerinfo.vcpu,
                    'MAX_CONSIST': self.MAX_CONSIST,
                    'randomize': self.randomize,
                    'nbases': self.n_bases,
                    }
            )
            # Reverse:
            self.ex(
                command=command_R,
                input_targets=input_targets_R,
                output_targets=output_targets_R,
                extra_params={
                    'vcpu': self.containerinfo.vcpu,
                    'MAX_CONSIST': self.MAX_CONSIST,
                    'randomize': self.randomize,
                    'nbases': self.n_bases,
                    }
            )


class DADA2_DADA(sl.ContainerTask):
    # DADA2 filter and trim for a specimen (F and R)
    container = 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'

    # Dependencies
    in_derep = None
    in_errM = None

    # Parameters
    specimen = sl.Parameter()
    path = sl.Parameter()

    def out_rds(self):
        rds_dict = {}
        rds_dict['R1'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R1.dada.rds".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        rds_dict['R2'] = sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.R2.dada.rds".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
        return rds_dict
    
    def run(self):
        input_targets = {
            'derep_1': self.in_derep()['R1'],
            'derep_2': self.in_derep()['R2'],
            'errM_1': self.in_errM()['R1'],
            'errM_2': self.in_errM()['R2'],
        }
        output_targets = {
            'dada_1': self.out_rds()['R1'],
            'dada_2': self.out_rds()['R2']
        }

        self.ex(
            command=(
                'Rscript -e "'
                "library('dada2'); "
                "errM_1 <- readRDS('$errM_1'); "
                "derep_1 <- readRDS('$derep_1'); "
                "dadaResult_1 <- dada(derep_1, err=errM_1, multithread=$vcpu, verbose=FALSE); "
                "saveRDS(dadaResult_1, '$dada_1'); "
                "errM_2 <- readRDS('$errM_2'); "
                "derep_2 <- readRDS('$derep_2'); "
                "dadaResult_2 <- dada(derep_2, err=errM_2, multithread=$vcpu, verbose=FALSE); "
                "saveRDS(dadaResult_2, '$dada_2'); "
                '"'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={'vcpu': self.containerinfo.vcpu}
        )


class DADA2_Merge(sl.ContainerTask):
    # DADA2 Merge forward and reverse reads
    container = 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'

    # Dependencies
    in_derep = None
    in_dada = None

    # Parameters
    specimen = sl.Parameter()
    path = sl.Parameter()

    def out_rds(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.merge.rds".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
    
    def run(self):
        input_targets = {
            'derep_1': self.in_derep()['R1'],
            'derep_2': self.in_derep()['R2'],
            'dada_1': self.in_dada()['R1'],
            'dada_2': self.in_dada()['R2'],
        }

        output_targets = {
            'merger': self.out_rds(),
        }

        self.ex(
            command=(
                'Rscript -e "'
                "library('dada2'); "
                "dada_1 <- readRDS('$dada_1'); "
                "derep_1 <- readRDS('$derep_1'); "
                "dada_2 <- readRDS('$dada_2'); "
                "derep_2 <- readRDS('$derep_2'); "
                'merger <- mergePairs('
                'dada_1, derep_1, '
                'dada_2, derep_2 '
                ', verbose=FALSE); '
                "saveRDS(merger, '$merger'); "
                '"'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
        )

class DADA2_Specimen_Seqtab(sl.ContainerTask):
    # DADA2 Make a sequence table for a specimen
    container = 'golob/dada2:1.6.0__bcw.0.3.0'

    # Dependencies
    in_merge = None

    # Parameters
    specimen = sl.Parameter()
    path = sl.Parameter()

    def out_rds(self):
        return sl.ContainerTargetInfo(
            self,
            os.path.join(
                self.path,
                "{}.seqtab.rds".format(self.specimen)
            ),
            format=luigi.format.Nop
        )
    
    def run(self):
        input_targets = {
            'merged': self.in_merge(),
        }

        output_targets = {
            'seqtab': self.out_rds(),
        }

        self.ex(
            command=(
                'Rscript -e "'
                "library('dada2'); "
                "merged <- readRDS('$merged'); "
                'seqtab <- makeSequenceTable(merged); '
                "rownames(seqtab) <- c('$specimen'); "  # Inject specimen name to ease later combination
                "saveRDS(seqtab, '$seqtab'); "
                '"'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'specimen': self.specimen,
            }
        )


class DADA2_Combine_Seqtabs_Lin(sl.ContainerTask):
    # DADA2 Combine seqtabs
    container = 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'

    # Dependencies
    in_seqtabs = None

    # Parameters
    fn = sl.Parameter()

    def out_rds(self):
        return sl.ContainerTargetInfo(
            self,
            self.fn,
            format=luigi.format.Nop
        )
    
    def run(self):
        # Strategy here is to loop through the seqtabs and slowly build up
        # the combined table

        # Initialize
        with self.out_rds().open('w') as out_h:
            out_h.write(b"")

        # Output target never changes
        input_targets = {
            'combined_seqtab': self.out_rds(),
        }

        for seq_tab_T in self.in_seqtabs:
            input_targets['indv_seqtab'] =  seq_tab_T()
            command = (
                'Rscript -e "'
                "library('dada2'); "
                "indv_seqtab <- readRDS('$indv_seqtab'); "
                'prior_seqtab <- tryCatch('
                    "{ readRDS('$combined_seqtab'); }, "
                    'error=function(cond){ print(cond); return (NA); }'
                '); '
                'valid_indv <- ('
                'is.matrix(indv_seqtab) && '
                'all(indv_seqtab)>=0 && '
                '!is.null(colnames(indv_seqtab)) && '
                '!is.null(rownames(indv_seqtab)) && '
                'all(sapply(colnames(indv_seqtab), nchar)>0) && '
                'all(sapply(rownames(indv_seqtab), nchar)>0) '
                '); '
                'valid_prior <- ('
                '!is.na(prior_seqtab) && '
                'is.matrix(prior_seqtab) && '
                'all(prior_seqtab)>=0 && '
                '!is.null(colnames(prior_seqtab)) && '
                '!is.null(rownames(prior_seqtab)) && '
                'all(sapply(colnames(prior_seqtab), nchar)>0) && '
                'all(sapply(rownames(prior_seqtab), nchar)>0) '
                '); '
                ' if (valid_prior && valid_indv) {'
                'combined_seqtab <- mergeSequenceTables(prior_seqtab, indv_seqtab); '
                "saveRDS(combined_seqtab, '$combined_seqtab') "
                '}'
                ' else if (valid_indv) { '
                "saveRDS(indv_seqtab, '$combined_seqtab') "
                '}'
                ' else if (valid_prior) { '
                "saveRDS(prior_seqtab, '$combined_seqtab') "
                '}'
                '"'
            )
            self.ex(
                command=command,
                input_targets=input_targets,
                inputs_mode='rw'
                #output_targets=output_targets,
            )
class DADA2_Combine_Seqtabs(sl.ContainerTask):
    # DADA2 Combine seqtabs
    container = 'golob/dada2-fast-combineseqtab:0.2.0_BCW_0.30A'

    # Dependencies
    in_seqtabs = None

    # Parameters
    fn = sl.Parameter()

    def out_rds(self):
        return sl.ContainerTargetInfo(
            self,
            self.fn,
            format=luigi.format.Nop
        )
    
    def run(self):
        input_targets = {"seqtab_{}".format(s_i): seqtab_t()
            for s_i, seqtab_t in enumerate(self.in_seqtabs)
        }

        output_targets = {
            'combined_seqtab': self.out_rds(),
        }

        command = (
            'combine_seqtab '
            '--rds $combined_seqtab '
            '--seqtabs '
        )

        command += " ".join(
            [
                "'${}'".format(st_id) 
            
            for st_id in input_targets.keys()]
        )

        self.ex(
            command=command,
            input_targets=input_targets,
            output_targets=output_targets,
        )

class DADA2_Combine_Seqtabs_Native(sl.ContainerTask):
    # DADA2 Combine seqtabs
    container = 'golob/dada2:1.8.0__bcw.0.3.0'

    # Dependencies
    in_seqtabs = None

    # Parameters
    fn = sl.Parameter()

    def out_rds(self):
        return sl.ContainerTargetInfo(
            self,
            self.fn,
            format=luigi.format.Nop
        )
    
    def run(self):
        input_targets = {"seqtab_{}".format(s_i): seqtab_t()
            for s_i, seqtab_t in enumerate(self.in_seqtabs)
        }

        output_targets = {
            'combined_seqtab': self.out_rds(),
        }

        command = (
            'Rscript -e "'
            "library('dada2'); "
            'seqtab_fns <- c('
        )

        command += ", ".join(
            [
                "'${}'".format(st_id) 
            
            for st_id in input_targets.keys()]
        )
        command += '); '
        # Apply to each element in the list a load operation. Store into seqtabs list
        command += 'seqtabs <- lapply(seqtab_fns, readRDS); '
        # Filter out tables that are empty (have no sequences), etc
        command += (
            'seqtabs.filtered <- '
            'seqtabs[lapply(seqtabs, '
            'function(tab) is.matrix(tab) '
            '&& all(tab>=0) '
            '&& !is.null(colnames(tab)) '
            '&& !is.null(rownames(tab)) '
            '&& all(sapply(colnames(tab), nchar)>0) '
            '&& all(sapply(rownames(tab), nchar)>0) ) == TRUE]; '
        )
        # Use the dada2 method to actually combine together
        command += 'combined_seqtab <- do.call(mergeSequenceTables, seqtabs.filtered); ' 
        command += "saveRDS(combined_seqtab, '$combined_seqtab') " # Save to disk
        # Closing quote for Rscript
        command += '"'

        self.ex(
            command=command,
            input_targets=input_targets,
            output_targets=output_targets,
        )


class DADA2_Remove_Chimera(sl.ContainerTask):
    # DADA2 Remove chimeras from Seqtab
    container = 'golob/dada2:1.8.0.ub.1804__bcw.0.3.0A'

    # Dependencies
    in_seqtab = None

    # Parameters
    fn_rds = sl.Parameter()
    fn_csv = sl.Parameter()
    method = sl.Parameter(default='consensus')

    def out_rds(self):
        return sl.ContainerTargetInfo(
            self,
            self.fn_rds,
            format=luigi.format.Nop
        )

    def out_csv(self):
        return sl.ContainerTargetInfo(
            self,
            self.fn_csv,
        )

    def run(self):
        input_targets = {
            'seqtab': self.in_seqtab()
        }

        output_targets = {
            'seqtab_nochim_rds': self.out_rds(),
            'seqtab_nochim_csv': self.out_csv(),
        }

        self.ex(
            command=(
                'Rscript -e "'
                "library('dada2'); "
                "seqtab <- readRDS('$seqtab'); "
                'seqtab_nochim <- removeBimeraDenovo('
                'seqtab, '
                "method='$method', "
                'multithread=$vcpu); '
                'print((sum(seqtab) - sum(seqtab_nochim)) / sum(seqtab)); '
                "saveRDS(seqtab_nochim, '$seqtab_nochim_rds'); "
                "write.csv(seqtab_nochim, '$seqtab_nochim_csv', na='') "
                '"'
            ),
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'method': self.method,
                'vcpu': self.containerinfo.vcpu
            }
        )


class DADA2_SV_to_PPlacer(sl.ContainerTask):
    # Convert DADA2 sequence variant table to outputs compatible with pplacer
    # guppy and rppr.
    container = 'golob/dada2-pplacer:0.2.0__bcw__0.3.0'

    # Dependencies
    in_seqtab_csv = None

    # Parameters
    fasta_fn = sl.Parameter()
    weights_fn = sl.Parameter()
    map_fn = sl.Parameter()

    def out_fasta(self):
        return sl.ContainerTargetInfo(
            self,
            self.fasta_fn,
        )

    def out_weights(self):
        return sl.ContainerTargetInfo(
            self,
            self.weights_fn,
        )

    def out_map(self):
        return sl.ContainerTargetInfo(
            self,
            self.map_fn,
        )

    def run(self):
        input_targets = {
            'seqtab': self.in_seqtab_csv()
        }

        output_targets = {
            'fasta': self.out_fasta(),
            'weights': self.out_weights(),
            'map': self.out_map(),
        }

        self.ex(
            command=(
                'dada2-seqtab-to-pplacer '
                '--seqtable $seqtab '
                '--fasta_out_sequences $fasta '
                '--weights $weights '
                '--map $map'

            ),
            input_targets=input_targets,
            output_targets=output_targets,
        )
