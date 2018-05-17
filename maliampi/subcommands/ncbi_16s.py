#!/usr/bin/python
import luigi
import sciluigi as sl
from lib.tasks import NT_AccessionsForQuery, NT_Repo_Update_Accessions, LoadFile, NT_Repo_Fill

import os

ENGINE = 'docker'


# Workflow
class Workflow_NCBI_16s(sl.WorkflowTask):
    #
    # Take a set of sequence variants in FASTA format and at least one repository
    # of reference sequences.
    # Search the repository / repositories for matches above a specified threshold
    # for the sequence variants.
    #  Use those recruited full length repo sequences to build a refpkg.
    #
    working_dir = sl.Parameter()
    ncbi_email = sl.Parameter()
    repo_url = sl.Parameter()

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
        # Load current accessions with 16s in a genome
        #

        repo_url = self.new_task(
            'load_repo_url',
            LoadFile,
            path=self.repo_url,
        )

        acc_genome_16s = self.new_task(
            'genome_16s_accessions',
            NT_AccessionsForQuery,
            containerinfo=self.test_containerinfo,
            email=self.ncbi_email,
            accessions_fn=os.path.join(
                self.working_dir,
                'ncbi_16s',
                'accession',
                'genome_16s.csv'
            ),
            query=(
                "16s[All Fields] AND rRNA[Feature Key]"
                " AND Bacteria[Organism]"
                " AND 500000 : 99999999999[Sequence Length]"
                " AND genome[All Fields]"
            ),
        )

        repo_genome_update = self.new_task(
            'repo_genome_update',
            NT_Repo_Update_Accessions,
            extra_values={'is_genome': True},
        )
        repo_genome_update.in_repo_url = repo_url.out_file
        repo_genome_update.in_accessions = acc_genome_16s.out_accessions

        repo_filled = self.new_task(
            'repo_fill',
            NT_Repo_Fill,
            containerinfo=self.test_containerinfo,
            email=self.ncbi_email,
            working_dir=os.path.join(
                self.working_dir,
                'ncbi_16s',
            ),
        )
        repo_filled.in_repo = repo_genome_update.out_repo

        return(repo_filled, repo_genome_update)


def build_args(parser):
    parser.add_argument(
        '--ncbi-email',
        help="""email to use with NCBI""",
        required=True
        )
    parser.add_argument(
        '--repo-secret',
        help="""Path to a file containing the URL / secret for the NCBI repo
        """,
        type=str,
        required=True,
    )
    parser.add_argument(
        '--working-dir',
        help="""Path of a suitable working directory
        (defaults to the current working directory)""",
        type=str,
        default='.',
    )
