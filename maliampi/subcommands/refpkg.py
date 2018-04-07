#!/usr/bin/python
import luigi
import sciluigi as sl
from lib.tasks import LoadFastaSeqs, SearchRepoForMatches, ListVsearchOptions
from lib.containertask import ContainerInfo
import os


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
        # Search the sequence variants in filtered repository
        #

        search_sv_filtered = self.new_task(
            'search_sv_filtered',
            SearchRepoForMatches,
            matches_uc_path=os.path.join(self.working_dir,
                                         'repo_matches.types.uc'),
            unmatched_exp_seqs_path=os.path.join(self.working_dir, 
                                                 'exp_seqs_unmatched.types.fasta'),
            matched_repo_seqs_path=os.path.join(self.working_dir, 
                                                'recruited_repo_seqs.types.fasta'),
            min_id=self.min_id_types,
            maxaccepts=10,  # likewise, should be a parameter in a config file. For, take the top 10 (roughly corresponding to a 95% id for most)
        )
        search_sv_filtered.in_exp_seqs = sequence_variants.out_seqs
        search_sv_filtered.in_repo_seqs = repo_seqs_filtered.out_seqs

        return(search_sv_filtered)


def build_args(parser):
    parser.add_argument(
        '--sequence-variants',
        help="""Path to sequence variants (in FASTA format)
            for which we need a reference set created""",
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
