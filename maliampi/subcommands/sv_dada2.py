#!/usr/bin/python
import luigi
import sciluigi as sl
from lib.tasks import LoadManifest, LoadSpecimenReads, BCCSpecimenReads
from lib.tasks import DADA2_FilterAndTrim, DADA2_Dereplicate, DADA2_LearnError
from lib.tasks import DADA2_DADA, DADA2_Merge, DADA2_Specimen_Seqtab
from lib.tasks import DADA2_Combine_Seqtabs, DADA2_Remove_Chimera, DADA2_SV_to_PPlacer

from collections import defaultdict
import logging
import os

ENGINE = 'singularity_pbs'
log = logging.getLogger('sciluigi-interface')


# Workflow
class Workflow_DADA2(sl.WorkflowTask):
    #
    #  Take a suitable reference package and a set of sequence variants
    #  Place onto the maximum likelihood tree and return a jplace-format
    #  along with some QC data.
    #  For now based on PPLACER, but with an option for others in the future
    #
    working_dir = sl.Parameter()
    destination_dir = sl.Parameter()
    manifest = sl.Parameter()
    barcodecop = sl.Parameter(default=True)

    test_containerinfo = sl.ContainerInfo(
                vcpu=2,
                mem=4096,
                container_cache='/scratch/schmidti_fluxm/golobj/singularity_img/',
                engine=ENGINE,
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
                slurm_partition='boneyard',
                pbs_account='schmidti_fluxm',
                pbs_queue='fluxm'
            )

    local_containerinfo = sl.ContainerInfo(
                vcpu=1,
                mem=4096,
                container_cache='/scratch/schmidti_fluxm/golobj/singularity_img/',
                engine=ENGINE,
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
                slurm_partition='boneyard',
                pbs_account='schmidti_fluxm',
                pbs_queue='fluxm',
                pbs_scriptpath='/scratch/schmidti_fluxm/golobj/scripttmp',
                timeout=10,
            )

    light_containerinfo = sl.ContainerInfo(
                vcpu=1,
                mem=4096,
<<<<<<< HEAD
                container_cache=os.path.abspath(os.path.join('../working', 'containers/')),
                #engine='aws_batch',
                engine='docker',
=======
                container_cache='/scratch/schmidti_fluxm/golobj/singularity_img/',
                engine=ENGINE,
>>>>>>> ae7dc910da3b93b93ca387f4ac951fa03e564177
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
                slurm_partition='boneyard',
                pbs_account='schmidti_fluxm',
                pbs_queue='fluxm',
                pbs_scriptpath='/scratch/schmidti_fluxm/golobj/scripttmp',
                timeout=10,
            )

    long_containerinfo = sl.ContainerInfo(
                vcpu=1,
                mem=4096,
                container_cache='/scratch/schmidti_fluxm/golobj/singularity_img/',
                engine=ENGINE,
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
                slurm_partition='boneyard',
                pbs_account='schmidti_fluxm',
                pbs_queue='fluxm',
                pbs_scriptpath='/scratch/schmidti_fluxm/golobj/scripttmp',
                timeout=6000,
            )

    heavy_containerinfo = sl.ContainerInfo(
                vcpu=6,
                mem=6 * 4096,
                container_cache='/scratch/schmidti_fluxm/golobj/singularity_img/',
                engine=ENGINE,
                aws_s3_scratch_loc='s3://fh-pi-fredricks-d/lab/golob/sl_temp/',
                aws_jobRoleArn='arn:aws:iam::064561331775:role/fh-pi-fredricks-d-batchtask',
                aws_batch_job_queue='optimal',
                slurm_partition='boneyard',
                pbs_account='schmidti_fluxm',
                pbs_queue='fluxm',
                pbs_scriptpath='/scratch/schmidti_fluxm/golobj/scripttmp',
                timeout=60,
            )

    heavy_containerinfo = light_containerinfo
    
    def workflow(self):
        #
        #  Load the manifest of files
        #
        manifest = self.new_task(
            'load_manifest',
            LoadManifest,
            path=self.manifest,
        )

        # For each specimen....
        specimen_tasks = defaultdict(dict)
        specimens = manifest.get_specimens()
        for specimen in specimens:
            # Load the specimen reads.
            specimen_tasks[specimen]['reads'] = self.new_task(
                'specimen_load_{}'.format(specimen),
                LoadSpecimenReads,
                specimen=specimen
            )
            specimen_tasks[specimen]['reads'].in_manifest = manifest.out_file
            if self.barcodecop and manifest.has_index() and manifest.is_paired():
                specimen_tasks[specimen]['verified_reads'] = self.new_task(
                    'specimen_bcc_{}'.format(specimen),
                    BCCSpecimenReads,
                    containerinfo=self.light_containerinfo,
                    specimen=specimen,
                    path=os.path.join(
                        self.working_dir,
                        'sv',
                        'bcc'
                    )
                )
                specimen_tasks[specimen]['verified_reads'].in_reads = specimen_tasks[specimen]['reads'].out_reads
            else:
                specimen_tasks[specimen]['verified_reads'] = specimen_tasks[specimen]['reads']

            # DADA2 filer and trim
            specimen_tasks[specimen]['dada2_ft'] = self.new_task(
                'dada2_ft_{}'.format(specimen),
                DADA2_FilterAndTrim,
                containerinfo=self.light_containerinfo,
                specimen=specimen,
                f_trunc=235,
                r_trunc=235,
                trim_left=15,
                path=os.path.join(
                    self.working_dir,
                    'sv',
                    'dada2',
                    'ft'
                )
            )
            specimen_tasks[specimen]['dada2_ft'].in_reads = specimen_tasks[specimen]['verified_reads'].out_reads

            specimen_tasks[specimen]['dada2_derep'] = self.new_task(
                'dada2_derep_{}'.format(specimen),
                DADA2_Dereplicate,
                containerinfo=self.light_containerinfo,
                specimen=specimen,
                path=os.path.join(
                    self.working_dir,
                    'sv',
                    'dada2',
                    'derep'
                )
            )
            specimen_tasks[specimen]['dada2_derep'].in_reads = specimen_tasks[specimen]['dada2_ft'].out_reads

        # Now we need the specimens grouped by batch to create error models. 
        batch_errModels = {}
        for batch, batched_specimens in manifest.batched_specimens():
            batch_errModels[batch] = self.new_task(
                'dada2_learn_error_batch_{}'.format(batch),
                DADA2_LearnError,
                containerinfo=self.heavy_containerinfo,
                batch=batch,
                tar_reads=False,
                path=os.path.join(
                    self.working_dir,
                    'sv',
                    'dada2',
                    'errM'
                )
            )
            batch_errModels[batch].in_reads = [
                specimen_tasks[s]['dada2_ft'].out_reads
                for s in specimen_tasks
                if s in batched_specimens
            ]
            for specimen in batched_specimens:
                specimen_tasks[specimen]['dada2_errM'] = batch_errModels[batch]

        # Back to for each specimen...
        for specimen in specimens:
            # DADA
            specimen_tasks[specimen]['dada2_dada'] = self.new_task(
                'dada2_dada_{}'.format(specimen),
                DADA2_DADA,
                containerinfo=self.heavy_containerinfo,
                specimen=specimen,
                path=os.path.join(
                    self.working_dir,
                    'sv',
                    'dada2',
                    'dada'
                )
            )
            specimen_tasks[specimen]['dada2_dada'].in_derep = specimen_tasks[specimen]['dada2_derep'].out_rds
            specimen_tasks[specimen]['dada2_dada'].in_errM = specimen_tasks[specimen]['dada2_errM'].out_rds

            # MERGE
            specimen_tasks[specimen]['dada2_merge'] = self.new_task(
                'dada2_merge_{}'.format(specimen),
                DADA2_Merge,
                containerinfo=self.light_containerinfo,
                specimen=specimen,
                path=os.path.join(
                    self.working_dir,
                    'sv',
                    'dada2',
                    'merged'
                )
            )
            specimen_tasks[specimen]['dada2_merge'].in_dada = specimen_tasks[specimen]['dada2_dada'].out_rds
            specimen_tasks[specimen]['dada2_merge'].in_derep = specimen_tasks[specimen]['dada2_derep'].out_rds

            # Seqtab
            specimen_tasks[specimen]['dada2_seqtab'] = self.new_task(
                'dada2_seqtab_{}'.format(specimen),
                DADA2_Specimen_Seqtab,
                containerinfo=self.light_containerinfo,
                specimen=specimen,
                path=os.path.join(
                    self.working_dir,
                    'sv',
                    'dada2',
                    'seqtab'
                )
            )
            specimen_tasks[specimen]['dada2_seqtab'].in_merge = specimen_tasks[specimen]['dada2_merge'].out_rds

        # Combine seqtabs
        combined_seqtab = self.new_task(
            'dada2_combine_seqtabs',
            DADA2_Combine_Seqtabs,
            containerinfo=self.long_containerinfo,
            fn=os.path.join(
                        self.working_dir,
                        'sv',
                        'dada2',
                        'seqtab.combined.rds'
                    )
        )
        combined_seqtab.in_seqtabs = [
                specimen_tasks[s]['dada2_seqtab'].out_rds
                for s in specimens
            ]

        combined_seqtab_nochim = self.new_task(
            'dada2_remove_chimera',
            DADA2_Remove_Chimera,
            containerinfo=self.heavy_containerinfo,
            fn_rds=os.path.join(
                        self.working_dir,
                        'sv',
                        'dada2',
                        'seqtab.combined.nochim.rds'
                    ),
            fn_csv=os.path.join(
                        self.destination_dir,
                        'seqtab.combined.nochim.csv'
                    )
        )
        combined_seqtab_nochim.in_seqtab = combined_seqtab.out_rds

        dada2_sv_to_pplacer = self.new_task(
            'dada2_sv_to_pplacer',
            DADA2_SV_to_PPlacer,
            containerinfo=self.light_containerinfo,
            fasta_fn=os.path.join(
                        self.destination_dir,
                        'dada2.sv.fasta',
                    ),
            weights_fn=os.path.join(
                        self.destination_dir,
                        'dada2.sv.weights.csv',
                    ),
            map_fn=os.path.join(
                        self.destination_dir,
                        'dada2.sv.map.csv',
                    )
        )
        dada2_sv_to_pplacer.in_seqtab_csv = combined_seqtab_nochim.out_csv

        return (dada2_sv_to_pplacer)


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
        '-M', '--manifest',
        help="""Manifest of files in CSV format
        one row per specimen. Must at least have columns labeled specimen and read__1
        """,
        type=str,
        required=True,
    )

