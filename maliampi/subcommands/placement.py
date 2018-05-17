#!/usr/bin/python
import luigi
import sciluigi as sl
from lib.tasks import LoadFile, CMAlignSeqs, LoadFastaSeqs, AlignmentStoToFasta
from lib.tasks import ExtractRefpkgAlignment, CombineAlignmentsSTO
from lib.tasks import PPLACER_PlaceAlignment, Jplace_Reduplicate
from lib.tasks import Jplace_PCA, Jplace_ADCL, Jplace_EDPL, Jplace_KR_Distance
from lib.tasks import Jplace_Alpha_Diversity, LoadRefpkgTGZ

import os

ENGINE = 'docker'


# Workflow
class Workflow_Placement(sl.WorkflowTask):
    #
    #  Take a suitable reference package and a set of sequence variants
    #  Place onto the maximum likelihood tree and return a jplace-format
    #  along with some QC data.
    #  For now based on PPLACER, but with an option for others in the future
    #
    working_dir = sl.Parameter()
    destination_dir = sl.Parameter()
    sv_fasta = sl.Parameter()
    refpkg_tgz = sl.Parameter()
    seq_map_csv = sl.Parameter()
    sv_weights_csv = sl.Parameter(default=None)

    test_containerinfo = sl.ContainerInfo(
                vcpu=2,
                mem=4096,
                container_cache=os.path.abspath(os.path.join('../working', 'containers/')),
                engine=ENGINE,
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
                slurm_partition='boneyard'
            )

    def workflow(self):
        #
        #  Load the refpkg (in tgz format)
        #
        refpkg_tgz = self.new_task(
            'load_refpkg_tgz',
            LoadRefpkgTGZ,
            path=self.refpkg_tgz,
            file_format='gzip',
        )

        # Load the seq map
        seq_map = self.new_task(
            'load_seq_map',
            LoadFile,
            path=self.seq_map_csv
        )

        # Load the weights if provided
        if self.sv_weights_csv:
            sv_weights = self.new_task(
                'load_sv_weight',
                LoadFile,
                path=self.sv_weights_csv
            )
        else:
            sv_weights = None

        #  And unpack the refpkg to the relevant bits
        refpkg_alignments = self.new_task(
            'refpkg_alignments',
            ExtractRefpkgAlignment,
            aln_fasta_fn=os.path.join(
                self.working_dir,
                'placement',
                'refpkg.aln.fasta'
            ),
            aln_sto_fn=os.path.join(
                self.working_dir,
                'placement',
                'refpkg.aln.sto'
            ),
        )
        refpkg_alignments.in_refpkg_tgz = refpkg_tgz.out_refpkg_tgz

        #
        #  Load the sequence variants (fasta format)
        #
        sv_fasta = self.new_task(
            'load_sv',
            LoadFastaSeqs,
            fasta_seq_path=self.sv_fasta
        )

        #
        #  Align the sequence variants
        #
        sv_aligned = self.new_task(
            'align_sv',
            CMAlignSeqs,
            containerinfo=self.test_containerinfo,
            alignment_sto_fn=os.path.join(
                self.working_dir,
                'placement',
                'sv.aln.sto'
            ),
            alignment_score_fn=os.path.join(
                self.working_dir,
                'placement',
                'sv.aln.scores'
            ),
        )
        sv_aligned.in_seqs = sv_fasta.out_seqs

        sv_aligned_fasta = self.new_task(
            'align_sv_to_fasta',
            AlignmentStoToFasta,
            align_fasta_fn=os.path.join(
                self.working_dir,
                'placement',
                'sv.aln.fasta'
            ),
        )
        sv_aligned_fasta.in_align_sto = sv_aligned.out_align_sto

        #
        #  Combine the refpkg alignment with the sequence variant alignment
        #

        sv_refpkg_aln_sto = self.new_task(
            'combine_sv_refpkg_aln_sto',
            CombineAlignmentsSTO,
            containerinfo=self.test_containerinfo,
            combined_aln_sto_fn=os.path.join(
                self.working_dir,
                'placement',
                'sv_refpkg_aln.sto'
            )
        )
        sv_refpkg_aln_sto.in_aln_sto_1 = refpkg_alignments.out_aln_sto
        sv_refpkg_aln_sto.in_aln_sto_2 = sv_aligned.out_align_sto

        #
        #  Place the sequence variants using this combined aligment
        #
        dedup_jplace = self.new_task(
            'make_dedup_jplace',
            PPLACER_PlaceAlignment,
            containerinfo=self.test_containerinfo,
            jplace_fn=os.path.join(
                self.destination_dir,
                'placement',
                'dedup.jplace'
            )
        )
        dedup_jplace.in_refpkg_tgz = refpkg_tgz.out_refpkg_tgz
        dedup_jplace.in_merged_aln_sto = sv_refpkg_aln_sto.out_aln_sto

        #
        #  Reduplicate
        #

        if not sv_weights:
            redup_jplace = dedup_jplace
        else:
            redup_jplace = self.new_task(
                'reduplicate_jplace',
                Jplace_Reduplicate,
                containerinfo=self.test_containerinfo,
                jplace_fn=os.path.join(
                    self.destination_dir,
                    'placement',
                    'redup.jplace.gz'
                )
            )
            redup_jplace.in_jplace = dedup_jplace.out_jplace
            redup_jplace.in_weights = sv_weights.out_file

        #
        #  ADCL
        #
        adcl = self.new_task(
            'create_adcl',
            Jplace_ADCL,
            containerinfo=self.test_containerinfo,
            adcl_fn=os.path.join(
                self.destination_dir,
                'placement',
                'adcl.gz'
            )
        )
        adcl.in_jplace = redup_jplace.out_jplace

        #
        #  EDPL
        #

        edpl = self.new_task(
            'calculate_edpl',
            Jplace_EDPL,
            containerinfo=self.test_containerinfo,
            edpl_fn=os.path.join(
                self.destination_dir,
                'placement',
                'edpl.gz'
            )
        )
        edpl.in_jplace = redup_jplace.out_jplace

        #
        #  EPCA
        #
        epca = self.new_task(
            'calculate_epca',
            Jplace_PCA,
            containerinfo=self.test_containerinfo,
            path=os.path.join(
                self.destination_dir,
                'placement',
                'pca'
            ),
            prefix='epca',
            pca='epca'
        )
        epca.in_refpkg_tgz = refpkg_tgz.out_refpkg_tgz
        epca.in_seq_map = seq_map.out_file
        epca.in_jplace = redup_jplace.out_jplace

        #
        #  LPCA
        #

        lpca = self.new_task(
            'calculate_lpca',
            Jplace_PCA,
            containerinfo=self.test_containerinfo,
            path=os.path.join(
                self.destination_dir,
                'placement',
                'pca'
            ),
            prefix='lpca',
            pca='lpca'
        )
        lpca.in_refpkg_tgz = refpkg_tgz.out_refpkg_tgz
        lpca.in_seq_map = seq_map.out_file
        lpca.in_jplace = redup_jplace.out_jplace

        #
        #  KR-distance
        #

        kr_distance = self.new_task(
            'calculate_kr_distance',
            Jplace_KR_Distance,
            containerinfo=self.test_containerinfo,
            kr_fn=os.path.join(
                self.destination_dir,
                'placement',
                'kr_distance.csv'
            ),
        )
        kr_distance.in_refpkg_tgz = refpkg_tgz.out_refpkg_tgz
        kr_distance.in_seq_map = seq_map.out_file
        kr_distance.in_jplace = redup_jplace.out_jplace

        # 
        #  Alpha-Diversity
        #

        alpha_diversity = self.new_task(
            'calculate_alpha_diversity',
            Jplace_Alpha_Diversity,
            containerinfo=self.test_containerinfo,
            alpha_diversity_fn=os.path.join(
                self.destination_dir,
                'placement',
                'alpha_diversity.csv'
            ),
        )
        alpha_diversity.in_refpkg_tgz = refpkg_tgz.out_refpkg_tgz
        alpha_diversity.in_seq_map = seq_map.out_file
        alpha_diversity.in_jplace = redup_jplace.out_jplace

        return(epca, lpca, adcl, edpl, kr_distance, alpha_diversity)


def build_args(parser):
    parser.add_argument(
        '--working-dir',
        help="""Path of a suitable working directory
        (defaults to the current working directory)""",
        type=str,
        default='.',
    )
    parser.add_argument(
        '--destination-dir',
        help="""Path of a suitable destination directory
        for the various placement outputs""",
        type=str,
        required=True,
    )
    parser.add_argument(
        '-sv', '--sequence-variants',
        help="""FASTA file to be placed sequence variants
        """,
        type=str,
        required=True,
    )
    parser.add_argument(
        '-rpkg', '--refpkg-tgz',
        help="""Tar-gzipped file containing the reference package
        to be placed against
        """,
        type=str,
        required=True
    )
    parser.add_argument(
        '-seq-map', '--seq-map-csv',
        help="""Map of sequence IDs to specimens, headerless csv seqid, specimen
        """,
        type=str,
    )
    parser.add_argument(
        '-weights', '--sv-weights-csv',
        help="""Weights of sequence variants for reduplication.
        Headerless csv: sv_id, seq_id, weight
        (optional)
        """,
        type=str,
    )
