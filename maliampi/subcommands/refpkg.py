#!/usr/bin/python
import luigi
import sciluigi as sl
from lib.tasks import LoadFastaSeqs, SearchRepoForMatches, CMAlignSeqs, RAxMLTree, LoadFile
from lib.tasks import BuildTaxtasticDB, FilterSeqinfoToFASTA, TaxTableForSeqInfo, ObtainCM
from lib.tasks import AlignmentStoToFasta, CombineRefpkg, CombineRepoMatches, ConfirmSeqInfoTaxonomy
from lib.tasks import CleanupTreeInfo
from lib.tasks import CombineRepoMatches
import os
import logging

log = logging.getLogger('sciluigi-interface')


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

    # Our goal is to have some annotated reference sequence with this
    # sequence identity for each experimental sequence variant.
    min_best = sl.Parameter()

    repo_seq_info = sl.Parameter()

    # Annotated sequences have trusted annotations
    # (taxonomic / gene content / etc). To be used for classification
    # this needs to be a subset where SOME metrics are used to validate
    # annotations (source, as in genome or type strain; consensus based on seq id)
    repo_annotated_fasta = sl.Parameter()

    min_id_annotated = sl.Parameter()

    # Not every experimental sequence variant is going to have an annotated sequence
    # availble. Thus we can have a larger second repo (where we do NOT trust the annotations)
    # but where at least the sequence quality is valid (full length, actually 16S, 
    # without many ambiguous bases) from which we can recruit additional sequences.
    # These can be annotated via some metric, OR used to to recruit additional seqs from the 
    # annotated set based on FULL LENGTH identity. 
    repo_valid_fasta = sl.Parameter(default="")
    min_id_valid = sl.Parameter()

    # email for entez
    entrez_email = sl.Parameter()

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
        log.info("Loaded sequence variants")

        # Load the sequence information
        seq_info_files = [
            self.new_task(
                'load_si_{}'.format(si_i),
                LoadFile,
                path=si_path
            )
            for si_i, si_path in enumerate(self.repo_seq_info.split(','))
        ]
        log.info("Loaded %d sequence information files", len(seq_info_files))

        #
        # Load the annotated repositories
        #
        repo_annotated = [
            self.new_task(
                'load_annotated_repo_{}'.format(r_i),
                LoadFastaSeqs,
                fasta_seq_path=r_path
            )
            for r_i, r_path in enumerate(self.repo_annotated_fasta.split(','))
        ]
        log.info("Loaded %d Annotated Repositories", len(repo_annotated))

        #
        # Search the sequence variants in the annotated repository
        #
        search_sv_annotated = []
        for ra_i, r_annotated in enumerate(repo_annotated):

            r_a_task = self.new_task(
                'search_sv_annotated_{}'.format(ra_i),
                SearchRepoForMatches,
                containerinfo=midcpu_containerinfo,
                matches_uc_path=os.path.join(self.working_dir,
                                            'refpkg',
                                            'repo.annotated__{}.matches.uc'.format(ra_i)),
                unmatched_exp_seqs_path=os.path.join(self.working_dir,
                                                    'refpkg',
                                                    'repo.annotated__{}.annotated.exp_seqs_unmatched.fasta'.format(ra_i)),
                matched_repo_seqs_path=os.path.join(self.working_dir,
                                                    'refpkg',
                                                    'repo.annotated__{}.recruited_repo_seqs.fasta'.format(ra_i)),
                min_id=self.min_id_annotated,
                maxaccepts=10,  # Default take the top 10 (roughly corresponding to a 95% id for most)
            )
            r_a_task.in_exp_seqs = sequence_variants.out_seqs
            r_a_task.in_repo_seqs = r_annotated.out_seqs
            search_sv_annotated.append(r_a_task)
        #
        # Combine Recruits into one file
        #

        combined_repo_matches = self.new_task(
            'combine_repo_matches',
            CombineRepoMatches,
            seqs_fn=os.path.join(self.working_dir,
                                        'refpkg',
                                        'combined.repo.maches.fasta'
            ),
            seq_info_fn=os.path.join(self.working_dir,
                                        'refpkg',
                                        'combined.repo.maches.seq_info.csv'
            ),
        )
        combined_repo_matches.in_seqs = [ssv.out_matched_repo_seqs for ssv in search_sv_annotated]
        combined_repo_matches.in_seq_info = [sif.out_file for sif in seq_info_files]

        refpkg_seqs = combined_repo_matches.out_seqs
        refpkg_seqinfo = combined_repo_matches.out_seq_info

        #
        # Verify the taxonomy for the refpkg seqinfo file.
        #

        verified_refpkg_seqinfo = self.new_task(
            'verify_refpkg_seqinfo_taxonomy',
            ConfirmSeqInfoTaxonomy,
            email=self.entrez_email,
            containerinfo=light_containerinfo,
            confirmed_seqinfo_path=os.path.join(
                self.working_dir,
                'refpkg',
                'seq_info.refpkg.verified_tax.csv'
            )
        )
        verified_refpkg_seqinfo.in_seq_info = refpkg_seqinfo
        verified_refpkg_seqinfo.in_tax_db = taxonomy_db.out_tax_db

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
        # Cleanup the tree info to remove cruft
        #

        tree_info_cleanup = self.new_task(
            'tree_info_cleanup',
            CleanupTreeInfo,
            tree_info_path=os.path.join(self.working_dir,
                                         'refpkg',
                                         'refpkg.tre.cleaned.info'),
        )
        tree_info_cleanup.in_tree_info = raxml_tree.out_tree_stats


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
        refpkg_taxtable.in_seq_info = verified_refpkg_seqinfo.out_seq_info
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
        combine_refpgk.in_tree_stats = tree_info_cleanup.out_tree_info
        combine_refpgk.in_taxtable = refpkg_taxtable.out_taxtable
        combine_refpgk.in_seq_info = verified_refpkg_seqinfo.out_seq_info
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
        '--entrez-email',
        help="Valid email for use with NCBI Entrez",
        required=True,
    )
    parser.add_argument(
        '--repo-seq-info',
        help="""Path to repository sequences information
            csv format expected. Multiple allowed.""",
        type=str,
        required=True,
        nargs='+'
    )
    parser.add_argument(
        '--repo-annotated-fasta',
        help="""Path(s) to repository sequences with trusted annotations
            from which we should recruit. FASTA format expected""",
        type=str,
        required=True,
        nargs='+',
    )
    parser.add_argument(
        '--repo-valid-fasta',
        help="""Path(s) to repository full-length sequences
            from which we should recruit after trying annotated first.
            FASTA format expected""",
        type=str,
        nargs='+',
        default=[],
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
        '--min-id-annotated',
        default=0.8,
        type=float,
        help='Min percent identity when recruiting from annotated-annotation 16S'
    )
    parser.add_argument(
        '--min-id-valid',
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
