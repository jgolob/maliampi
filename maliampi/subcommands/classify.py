#!/usr/bin/python
import luigi
import sciluigi as sl
from lib.tasks import LoadFile, CMAlignSeqs, LoadFastaSeqs, AlignmentStoToFasta
from lib.tasks import LoadRefpkgTGZ
from lib.tasks import ExtractRefpkgAlignment, CombineAlignmentsSTO
from lib.tasks import PlacementDB_Prep, PlacementDB_Classify_SV, PlacementDB_MCC
from lib.tasks import PlacementDB_AddSI
from lib.tasks import GenerateTables

import logging
import os

ENGINE = 'docker'
log = logging.getLogger('sciluigi-interface')

# Workflow
class Workflow_Classify(sl.WorkflowTask):
    #
    #  Take a suitable reference package and a set of sequence variants
    #  Place onto the maximum likelihood tree and return a jplace-format
    #  along with some QC data.
    #  For now based on PPLACER, but with an option for others in the future
    #
    working_dir = sl.Parameter()
    destination_dir = sl.Parameter()
    sv_fasta = sl.Parameter()
    jplace = sl.Parameter()

    refpkg_tgz = sl.Parameter()
    seq_map_csv = sl.Parameter()
    sv_weights_csv = sl.Parameter(default=None)
    labels = sl.Parameter(default=None)

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

        jplace = self.new_task(
            'load_jplace',
            LoadFile,
            path=self.jplace,
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

        if self.labels:
            labels = self.new_task(
                'load_labels',
                LoadFile,
                path=self.labels
            )
        else:
            labels = None

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
        #  Prep the placements.db using the refpkg
        #

        prepped_placementdb = self.new_task(
            'prep_placementdb',
            PlacementDB_Prep,
            containerinfo=self.test_containerinfo,
            placement_db_fn=os.path.join(
                self.destination_dir,
                'classification',
                'placement.db'
            )
        )
        prepped_placementdb.in_refpkg_tgz = refpkg_tgz.out_refpkg_tgz

        #
        #  Insert the seq_info / map of sv -> specimens
        #

        placement_db_w_si = self.new_task(
            'placement_db_add_si',
            PlacementDB_AddSI,
            containerinfo=self.test_containerinfo,
        )
        placement_db_w_si.in_placement_db = prepped_placementdb.out_placement_db
        placement_db_w_si.in_seq_map = seq_map.out_file

        #
        #  Classify the sequence variants
        #

        placement_db_classified = self.new_task(
            'classify_into_placement_db',
            PlacementDB_Classify_SV,
            containerinfo=self.test_containerinfo,
        )
        placement_db_classified.in_placement_db = placement_db_w_si.out_placement_db
        placement_db_classified.in_refpkg_tgz = refpkg_tgz.out_refpkg_tgz
        placement_db_classified.in_sv_refpkg_aln_sto = sv_refpkg_aln_sto.out_aln_sto
        placement_db_classified.in_jplace = jplace.out_file

        #
        #  Multiclass concat names
        #

        placement_db_mcc = self.new_task(
            'placement_db_multiclass_concat',
            PlacementDB_MCC,
            containerinfo=self.test_containerinfo,
        )
        placement_db_mcc.in_placement_db = placement_db_classified.out_placement_db
        placement_db_mcc.in_weights = sv_weights.out_file

        #
        #  Tabular CSV outputs
        #
        tables_for_rank = {}
        for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
            tables_for_rank[rank] = self.new_task(
                'by_specimen_{}'.format(rank),
                GenerateTables,
                containerinfo=self.test_containerinfo,
                tables_path=os.path.join(
                    self.destination_dir,
                    'classification',
                    'tables',
                ),
                rank=rank
            )
            tables_for_rank[rank].in_placement_db = placement_db_mcc.out_placement_db
            tables_for_rank[rank].in_seq_map = seq_map.out_file
            if labels:
                tables_for_rank[rank].in_labels = labels.out_file

        return (placement_db_mcc, tables_for_rank)


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
        '-JP', '--jplace',
        help="""Placements in jplace format
        """,
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
        required=True,
    )
    parser.add_argument(
        '-weights', '--sv-weights-csv',
        help="""Weights of sequence variants for reduplication.
        Headerless csv: sv_id, seq_id, weight
        (optional)
        """,
        type=str,
    )
    parser.add_argument(
        '-L', '--labels',
        help="""CSV of human-readable labels for specimens
        (optional)
        """,
        type=str,
    )
