import luigi
import sciluigi as sl
import os
from string import Template
from Bio import AlignIO, SeqIO
import shutil
from .targets import NCBI_Repo_Entries_TargetInfo,  NCBI_Repo_Filled_TargetInfo
from .targets import PlacementDB_Prepped_ContainerTargetInfo, RefpkgTGZ_ContainerTargetInfo
from .targets import PlacementDB_Classified_ContainerTargetInfo
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


class LoadManifest(LoadFile):
    
    def manifest(self, parameter_list):
        # Load manifest as csv and return as a list of dicts
        with self.out_file().open() as manifest_h:
            return [r for r in csv.DictReader(manifest_h)]


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
        from mem_top import mem_top
        for chunk_i, chunk_versions in enumerate(self.chunks(
                versions_need_fill,
                int(self.chunk_size))
                ):
            log.info("Updating chunk {} of {}".format(
                chunk_i+1,
                int(len(versions_need_fill) / self.chunk_size))
            )
            self.work_on_chunk(chunk_versions, gb_needed, raw_gb)
            log.info(mem_top())
            gc.collect()


class BuildTaxtasticDB(sl.ContainerTask):
    # A Task that uses vsearch to find matches for experimental sequences in a repo of sequences

    # Define the container (in docker-style repo format) to complete this task
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0'

    # Where to put the sqlite database
    tax_db_path = sl.Parameter()

    def out_tax_db(self):
        return sl.ContainerTargetInfo(self, self.tax_db_path)

    def run(self):
        # Get our host paths for inputs and outputs
        output_targets = {
            'tax_db': self.out_tax_db(),
        }

        self.ex(
            command='taxit ' +
                    ' new_database' +
                    ' $tax_db',
            output_targets=output_targets,
            )


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
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0'

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
    container = 'golob/pplacer:1.1alpha19rc_BCW_0.3.0'

    # Parameters to tell us where the refpkg should be placed,
    # as in: refpkg_path/refpkg_name.refpkg_version.refpkg
    refpkg_path = sl.Parameter()
    # what should we call this refpkg
    refpkg_name = sl.Parameter()
    # version, defaults to a timestamp of YYYYmmdd_HHMMSS
    refpkg_version = sl.Parameter(default=datetime.now().strftime('%Y%m%d_%H%M%S'))

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
            command='mkdir -p /scratch && cd /scratch && taxit' +
            ' create --locus 16S' +
            ' --package-name $refpkg_name --clobber' +
            ' --aln-fasta $aln_fasta ' +
            ' --aln-sto $aln_sto ' +
            ' --tree-file $tree ' +
            ' --tree-stats $tree_stats ' +
            ' --taxonomy $taxtable ' +
            ' --seq-info $seq_info ' +
            ' --profile $cm ' +
            ' && tar czvf $refpkg_tgz $refpkg_name/*',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={'refpkg_name': self.refpkg_name_ver()}
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
            contents = json.load(tar_h.extractfile(tar_contents_dict['CONTENTS.json']))
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

        self.ex(
            command='mkdir -p /refpkg && cd /refpkg && tar xzvf $refpkg_tgz -C /refpkg/ ' +
                    ' && pplacer -p --inform-prior --prior-lower $prior_lower --map-identity' +
                    ' -j $nproc' +
                    ' -c /refpkg/$refpkg_rel_path' +
                    ' $merged_aln' +
                    ' -o $jplace',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'prior_lower': self.pplacer_prior_lower,
                'nproc': self.containerinfo.vcpu,
                'refpkg_rel_path': refpkg_rel_path,
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

        self.ex(
            command='mkdir -p /refpkg && cd /refpkg && tar xzvf $refpkg_tgz -C /refpkg/ '
                    ' && mkdir -p /results '
                    ' && guppy $pca' +
                    ' $jplace:$seq_map ' +
                    ' -c /refpkg/$refpkg_rel_path ' + 
                    ' --out-dir /results --prefix $prefix' +
                    ' && cp /results/$prefix.xml $xml' + 
                    ' && cp /results/$prefix.proj $proj' + 
                    ' && cp /results/$prefix.trans $trans',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'refpkg_rel_path': refpkg_rel_path,
                'prefix': self.prefix,
                'pca': pca
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

        self.ex(
            command='mkdir -p /refpkg && cd /refpkg && tar xzvf $refpkg_tgz -C /refpkg/ '
                    ' && guppy kr' +
                    ' --list-out' +
                    ' -c /refpkg/$refpkg_rel_path ' +
                    ' $jplace:$seq_map ' +
                    ' -o $kr_distance',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'refpkg_rel_path': refpkg_rel_path,
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
        refpkg_rel_path = self.in_refpkg_tgz().get_refpkg_rel_path()
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
            extra_params={
                'refpkg_rel_path': refpkg_rel_path,
            }
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

        self.ex(
            command='mkdir -p /refpkg && cd /refpkg && tar xzvf $refpkg_tgz -C /refpkg/ '
                    ' && rppr prep_db' +
                    ' -c /refpkg/$refpkg_rel_path ' +
                    ' --sqlite $placement_db',
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'refpkg_rel_path': refpkg_rel_path,
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
            ' -c /refpkg/$refpkg_rel_path '
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
            }
        extra_params.update(self.guppy_parameter_to_dict())

        self.ex(
            command='mkdir -p /refpkg && cd /refpkg && tar xzvf $refpkg_tgz -C /refpkg/ &&' +
                    self.make_guppy_command(),
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
                'mkdir -p /working &&'
                ' python2 /usr/bin/classif_table.py'
                ' $placement_db'
                ' /working/by_taxon.csv'
                ' --rank $rank'
                ' --specimen-map $seq_map'
                ' --by-specimen /working/by_specimen.csv'
                ' --tallies-wide /working/tallies_wide.csv'
                ' && cp -r /working/by_taxon.csv $by_taxon'
                ' && cp -r /working/by_specimen.csv $by_specimen'
                ' && cp -r /working/tallies_wide.csv $tallies_wide'
        )

        if self.in_labels:
            command += ' --metadata-map $labels'

        self.ex(
            command=command,
            input_targets=input_targets,
            output_targets=output_targets,
            extra_params={
                'rank': self.rank,
            }
        )