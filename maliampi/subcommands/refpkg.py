#!/usr/bin/python
import luigi
import sciluigi as sl
from lib.tasks import LoadFastaSeqs, SearchRepoForMatches, CMAlignSeqs, RAxMLTree, LoadFile
from lib.tasks import BuildTaxtasticDB, FilterSeqinfoToFASTA, TaxTableForSeqInfo, ObtainCM
from lib.tasks import AlignmentStoToFasta, CombineRefpkg, CombineRepoMatches
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
    repo_fl_seq_info = sl.Parameter(default="")
    repo_fl_fasta = sl.Parameter(default="")
    repo_primary_seq_info = sl.Parameter()
    repo_primary_fasta = sl.Parameter()
    min_id_primary = sl.Parameter()
    min_id_fl = sl.Parameter()
    min_best = sl.Parameter()

    def workflow(self):
        # Intialize our container info
        light_containerinfo = sl.ContainerInfo()
        light_containerinfo.from_config(
            section='light'
        )
        midcpu_containerinfo = sl.ContainerInfo()
        midcpu_containerinfo.from_config(
            section='midcpu'
        )
        heavy_containerinfo = sl.ContainerInfo()
        heavy_containerinfo.from_config(
            section='heavy'
        )
        highmem_containerinfo = sl.ContainerInfo()
        highmem_containerinfo.from_config(
            section='highmem'
        )

        #
        # Build our taxonomy db
        #
        taxonomy_db = self.new_task(
            'taxonomy_db',
            BuildTaxtasticDB,
            containerinfo=light_containerinfo,
            tax_db_path=os.path.join(
                self.working_dir,
                'refpkg',
                'taxonomy.db'
                )
        )

        #
        # Load the sequence variants
        #
        sequence_variants = self.new_task(
            'load_sequence_variants',
            LoadFastaSeqs,
            fasta_seq_path=self.sequence_variants_path
        )

        #
        # Load the primary repository
        #
        repo_primary = self.new_task(
            'load_primary_repo',
            LoadFastaSeqs,
            fasta_seq_path=self.repo_primary_fasta
        )

        repo_primary_seq_info = self.new_task(
            'load_primary_seq_info',
            LoadFile,
            path=self.repo_primary_seq_info
        )


        #
        # Search the sequence variants in the primary repository
        #

        search_sv_primary = self.new_task(
            'search_sv_primary',
            SearchRepoForMatches,
            containerinfo=midcpu_containerinfo,
            matches_uc_path=os.path.join(self.working_dir,
                                         'refpkg',
                                         'repo_matches.primary.uc'),
            unmatched_exp_seqs_path=os.path.join(self.working_dir,
                                                 'refpkg',
                                                 'exp_seqs_unmatched.primary.fasta'),
            matched_repo_seqs_path=os.path.join(self.working_dir,
                                                'refpkg',
                                                'recruited_repo_seqs.primary.fasta'),
            min_id=self.min_id_primary,
            maxaccepts=10,  # Default take the top 10 (roughly corresponding to a 95% id for most)
        )
        search_sv_primary.in_exp_seqs = sequence_variants.out_seqs
        search_sv_primary.in_repo_seqs = repo_primary.out_seqs

        # 
        #  Filter the primary seq_info to be limited to entries for our recruits
        #

        filter_seqinfo_primary = self.new_task(
            'filter_si_primary',
            FilterSeqinfoToFASTA,
            filtered_seq_info_fn=os.path.join(
                self.working_dir,
                'refpkg',
                'repo_matches.primary.seq_info.csv'
            )
        )
        filter_seqinfo_primary.in_fasta = search_sv_primary.out_matched_repo_seqs
        filter_seqinfo_primary.in_seq_info = repo_primary_seq_info.out_file

        refpkg_seqs = search_sv_primary.out_matched_repo_seqs
        refpkg_seqinfo = filter_seqinfo_primary.out_seq_info

        #
        # Parse UC file to determine if we achieved our minimum-best goal
        # for each SV.
        #

        #
        # Align recruited repo seqs
        #

        align_recruits = self.new_task(
            'align_recruits',
            CMAlignSeqs,
            containerinfo=highmem_containerinfo,
            alignment_sto_fn=os.path.join(
                self.working_dir,
                'refpkg',
                'recruit.aln.sto'
            ),
            alignment_score_fn=os.path.join(
                self.working_dir,
                'refpkg',
                'recruit.aln.scores'
            ),
        )
        align_recruits.in_seqs = refpkg_seqs

        #
        # Make a fasta version of the alignment
        #

        align_fasta = self.new_task(
            'align_fasta',
            AlignmentStoToFasta,
            align_fasta_fn=os.path.join(
                self.working_dir,
                'refpkg',
                'recruit.aln.fasta'
            ),
        )
        align_fasta.in_align_sto = align_recruits.out_align_sto

        #
        #  Make a tree of the reference package sequences
        #

        raxml_tree = self.new_task(
            'raxml_tree',
            RAxMLTree,
            containerinfo=heavy_containerinfo,
            tree_path=os.path.join(self.working_dir,
                                   'refpkg',
                                   'refpkg.tre'),
            tree_stats_path=os.path.join(self.working_dir,
                                         'refpkg',
                                         'refpkg.tre.info'),
        )
        raxml_tree.in_align_fasta = align_fasta.out_align_fasta

        #
        #  Start to assemble the reference package at this point
        #

        # Taxtable
        refpkg_taxtable = self.new_task(
            'refpkg_taxtable',
            TaxTableForSeqInfo,
            containerinfo=light_containerinfo,
            taxtable_path=os.path.join(
                self.working_dir,
                'refpkg',
                'taxtable.csv'
            )
        )
        refpkg_taxtable.in_seq_info = refpkg_seqinfo
        refpkg_taxtable.in_tax_db = taxonomy_db.out_tax_db

        # Covariance Matrix
        obtain_cm = self.new_task(
            'obtain_cm',
            ObtainCM,
            containerinfo=light_containerinfo,
            cm_destination=os.path.join(
                self.working_dir,
                'refpkg',
                'rRNA_16S_SSU.cm'
            )
        )

        # And the actual combination step
        combine_refpgk = self.new_task(
            'combine_refpkg',
            CombineRefpkg,
            containerinfo=light_containerinfo,
            refpkg_path=os.path.join(
                self.new_refpkg_path,
                'refpkg',
            ),
            refpkg_name=self.new_refpkg_name,
        )
        combine_refpgk.in_aln_fasta = align_fasta.out_align_fasta
        combine_refpgk.in_aln_sto = align_recruits.out_align_sto
        combine_refpgk.in_tree = raxml_tree.out_tree
        combine_refpgk.in_tree_stats = raxml_tree.out_tree_stats
        combine_refpgk.in_taxtable = refpkg_taxtable.out_taxtable
        combine_refpgk.in_seq_info = refpkg_seqinfo
        combine_refpgk.in_cm = obtain_cm.out_cm

        return(combine_refpgk)



        #
        # Combine the sequences, avoiding duplicate sequences
        #

        combined_recruits = self.new_task(
            'combine_repo_recruits',
            CombineRepoMatches,
            seqs_fn=os.path.join(
                self.working_dir,
                'refpkg',
                'recruits.combined.fasta'
            ),
            seq_info_fn=os.path.join(
                self.working_dir,
                'refpkg',
                'recruits.combined.seq_info.csv'
            )
        )
        combined_recruits.in_seqs = [
            search_sv_genomes.out_matched_repo_seqs,
            search_sv_filtered.out_matched_repo_seqs,
        ]
        combined_recruits.in_seq_info = [
            repo_genomes_seq_info.out_file,
            repo_filtered_seq_info.out_file,
        ]






        refpkg_taxtable = self.new_task(
            'refpkg_taxtable',
            TaxTableForSeqInfo,
            containerinfo=self.local_containerinfo,
            taxtable_path=os.path.join(
                self.working_dir,
                'refpkg',
                'taxtable.csv'
            )
        )
        refpkg_taxtable.in_seq_info = combined_recruits.out_seq_info
        refpkg_taxtable.in_tax_db = taxonomy_db.out_tax_db

        obtain_cm = self.new_task(
            'obtain_cm',
            ObtainCM,
            containerinfo=self.local_containerinfo,
            cm_destination=os.path.join(
                self.working_dir,
                'refpkg',
                'rRNA_16S_SSU.cm'
            )
        )

        combine_refpgk = self.new_task(
            'combine_refpkg',
            CombineRefpkg,
            containerinfo=self.local_containerinfo,
            refpkg_path=os.path.join(
                self.working_dir,
                'refpkg',
            ),
            refpkg_name='test',
        )
        combine_refpgk.in_aln_fasta = align_fasta.out_align_fasta
        combine_refpgk.in_aln_sto = align_recruits.out_align_sto
        combine_refpgk.in_tree = raxml_tree.out_tree
        combine_refpgk.in_tree_stats = raxml_tree.out_tree_stats
        combine_refpgk.in_taxtable = refpkg_taxtable.out_taxtable
        combine_refpgk.in_seq_info = combined_recruits.out_seq_info
        combine_refpgk.in_cm = obtain_cm.out_cm

        return(combine_refpgk)


def build_args(parser):
    parser.add_argument(
        '--sequence-variants',
        help="""Path to sequence variants (in FASTA format)
            for which we need a reference set created""",
        required=True
        )
    parser.add_argument(
        '--repo-fl-seq-info',
        help="""Path to full-length repository sequences information
            csv format expected""",
        type=str,
        default=""
        )
    parser.add_argument(
        '--repo-fl-fasta',
        help="""Path(s) to repository full-length sequences
            from which we should recruit. FASTA format expected""",
        type=str,
        default=""
        )
    parser.add_argument(
        '--repo-primary-seq-info',
        help="""Path to repository sequences information
            csv format expected""",
        type=str,
        required=True
        )
    parser.add_argument(
        '--repo-primary-fasta',
        help="""Path(s) to repository sequences with primary annotations
            from which we should recruit. FASTA format expected""",
        type=str,
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
        '--min-id-primary',
        default=0.8,
        type=float,
        help='Min percent identity when recruiting from primary-annotation 16S'
    )
    parser.add_argument(
        '--min-id-fl',
        default=0.99,
        type=float,
        help='Min percent identity when recruiting from full-length 16S'
    )
    parser.add_argument(
        '--min-best',
        default=1.0,
        type=float,
        help='Min percent identity to consider a query matched'
    )
