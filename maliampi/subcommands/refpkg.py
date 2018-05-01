#!/usr/bin/python
import luigi
import sciluigi as sl
from lib.tasks import LoadFastaSeqs, SearchRepoForMatches, CMAlignSeqs, RAxMLTree, LoadFile
from lib.tasks import AlignmentStoToFasta
import os

ENGINE = 'docker'

# Workflow
class WorkflowMakeRefpkg(sl.WorkflowTask):
    #
    # Take a set of sequence variants in FASTA format and at least one repository
    # of reference sequences.
    # Search the repository / repositories for matches above a specified threshold
    # for the sequence variants.
    #  Use those recruited full length repo sequences to build a refpkg.
    #
    working_dir = sl.Parameter()
    sequence_variants_path = sl.Parameter()
    new_refpkg_path = sl.Parameter()
    new_refpkg_name = sl.Parameter()
    repo_seq_info = sl.Parameter()
    repo_seqs_filtered = sl.Parameter()
    min_id_types = sl.Parameter()
    min_id_filtered = sl.Parameter()
    min_id_unnamed = sl.Parameter()
    min_best = sl.Parameter()

    def workflow(self):
        #
        # Load the sequence variants
        #
        sequence_variants = self.new_task(
            'load_sequence_variants',
            LoadFastaSeqs,
            fasta_seq_path=self.sequence_variants_path
        )

        #
        # Load the filtered repository
        #
        repo_seqs_filtered = self.new_task(
            'load_repo_seqs_filtered',
            LoadFastaSeqs,
            fasta_seq_path=self.repo_seqs_filtered
        )

        # 
        # Load the repository seq_info (that has the taxonomic information, etc)
        #

        repo_seq_info = self.new_task(
            'load_repo_seq_info',
            LoadFile,
            path=self.repo_seq_info
        )

        #
        # Search the sequence variants in filtered repository
        #

        filtered_search_sv = self.new_task(
            'filtered_search_sv',
            SearchRepoForMatches,
            containerinfo=sl.ContainerInfo(
                vcpu=2,
                mem=4096,
                container_cache=os.path.join(self.working_dir, 'containers/'),
                engine=ENGINE,
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
            ),
            matches_uc_path=os.path.join(self.working_dir,
                                         'repo_matches.filtered.uc'),
            unmatched_exp_seqs_path=os.path.join(self.working_dir,
                                                 'exp_seqs_unmatched.filtered.fasta'),
            matched_repo_seqs_path=os.path.join(self.working_dir,
                                                'recruited_repo_seqs.filtered.fasta'),
            min_id=self.min_id_filtered,
            maxaccepts=10,  # Default take the top 10 (roughly corresponding to a 95% id for most)
        )
        filtered_search_sv.in_exp_seqs = sequence_variants.out_seqs
        filtered_search_sv.in_repo_seqs = repo_seqs_filtered.out_seqs

        #
        # Fill 'lonely' recruits (where only one species represented for a genus in the recruits
        #

        # fill_lonely_recruits = self.new_task(
        #    'fill_lonely_recruits',
        #    FillLonely,
        #)

        #
        # Align the recruited repo sequences
        #
        filtered_align_recruits = self.new_task(
            'filtered_align_recruits',
            CMAlignSeqs,
            containerinfo=sl.ContainerInfo(
                vcpu=2,
                mem=4096,
                container_cache=os.path.join(self.working_dir, 'containers/'),
                engine=ENGINE,
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
            ),
            alignment_sto_fn=os.path.join(
                self.working_dir,
                'filtered.aln.sto'
            ),
            alignment_score_fn=os.path.join(
                self.working_dir,
                'filtered.aln.scores'
            ),
        )
        filtered_align_recruits.in_seqs = filtered_search_sv.out_matched_repo_seqs

        #
        # Make a fasta version of the alignment
        #

        filtered_align_fasta = self.new_task(
            'filtered_align_fasta',
            AlignmentStoToFasta,
            align_fasta_fn=os.path.join(
                self.working_dir,
                'filtered.aln.fasta'
            ),
        )
        filtered_align_fasta.in_align_sto = filtered_align_recruits.out_align_sto

        raxml_tree = self.new_task(
            'raxml_tree',
            RAxMLTree,
            containerinfo=sl.ContainerInfo(
                vcpu=2,
                mem=4096,
                container_cache=os.path.join(self.working_dir, 'containers/'),
                engine=ENGINE,
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
            ),
            tree_path=os.path.join(self.working_dir,
                                   'refpkg.tre'),
            tree_stats_path=os.path.join(self.working_dir,
                                         'refpkg.tre.info'),
            raxml_working_dir=os.path.join(self.working_dir,
                                           'refpkg_workdir/'),
        )
        raxml_tree.in_align_fasta = filtered_align_fasta.out_align_fasta

        return(raxml_tree)


def build_args(parser):
    parser.add_argument(
        '--sequence-variants',
        help="""Path to sequence variants (in FASTA format)
            for which we need a reference set created""",
        required=True
        )
    parser.add_argument(
        '--repo-seq-info',
        help="""Path to repository sequences information
            csv format expected""",
        type=str,
        required=True
        )
    parser.add_argument(
        '--repo-seqs-filtered',
        help="""Path(s) to repository sequences with trusted annotations
            from which we should recruit. FASTA format expected""",
        type=str,
        #  nargs='+',
        required=True
        )
    parser.add_argument(
        '--refpkg-destdir',
        help='Directory where the new reference package should be placed',
        required=True
    )
    parser.add_argument(
        '--refpkg-name',
        help='Name of the new refpkg (must be a valid filesystem name)',
        required=True
    )
    parser.add_argument(
        '--working-dir',
        help="""Path of a suitable working directory
        (defaults to the current working directory)""",
        type=str,
        default='.',
    )
    parser.add_argument(
        '--min-id-types',
        default=0.8,
        type=float,
        help='Min percent identity when recruiting from type strain 16S'
    )
    parser.add_argument(
        '--min-id-filtered',
        default=0.9,
        type=float,
        help='Min percent identity when recruiting from filtered 16S'
    )
    parser.add_argument(
        '--min-id-unnamed',
        default=0.99,
        type=float,
        help='Min percent identity when recruiting from unnamed / no taxonomy 16S'
    )
    parser.add_argument(
        '--min-best',
        default=1.0,
        type=float,
        help='Min percent identity to consider a query matched'
    )
