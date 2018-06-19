#!/usr/bin/env python
import argparse
import logging
import sys
import luigi
import sciluigi as sl

from subcommands import refpkg, ncbi_16s, placement, classify, sv_dada2


class MALIAMPI:
    def __init__(self):
        parser = argparse.ArgumentParser(description="""
        Maximum Likelihood Amplicon Pipeline: An amplicon (PCR / 16S) microbiome pipeline.""")

        if len(sys.argv) < 2:
            parser.print_help()
        else:
            # Common options here
            parser.add_argument(
                '--luigi-manager',
                help="""Run with luigi's work manager daemon (default False)""",
                action='store_true',
            )

            parser.add_argument(
                '--workers',
                help="""How many concurrent workers to use""",
                type=int,
                default=1,
            )

            self.__subparsers__ = parser.add_subparsers(
                title='command',
                dest='command'
            )

            # repo_16s command
            subparser_ncbi_16s = self.__subparsers__.add_parser(
                'ncbi_16s',
                description="""Update a repository of 16S sequences from NCBI NT.
                """)
            ncbi_16s.build_args(subparser_ncbi_16s)

            # sequence_variants_dada2
            subparser_sv_dada2 = self.__subparsers__.add_parser(
                'sv_dada2',
                description="""Make sequence variants using DADA2
                """)
            sv_dada2.build_args(subparser_sv_dada2)

            # refpkg command
            subparser_refpkg = self.__subparsers__.add_parser(
                'refpkg',
                description="""Make a reference package
                                appropriate for pplacer or other pipelines.""")
            refpkg.build_args(subparser_refpkg)

            # placement command
            subparser_placement = self.__subparsers__.add_parser(
                'placement',
                description="""Place sequence variants on a reference package
                """)
            placement.build_args(subparser_placement)

            # classify command
            subparser_classify = self.__subparsers__.add_parser(
                'classify',
                description="""Classify sequence variants using a placement
                """)
            classify.build_args(subparser_classify)

            self.__args__ = parser.parse_args()
            if self.__args__.command == 'refpkg':
                self.refpkg()
            elif self.__args__.command == 'ncbi_16s':
                self.ncbi_16s()
            elif self.__args__.command == 'placement':
                self.placement()
            elif self.__args__.command == 'classify':
                self.classify()
            elif self.__args__.command == 'sv_dada2':
                self.sv_dada2()

    def sv_dada2(self):
        """Generate sequence variants using DADA2."""
        args = self.__args__

        if args.luigi_manager:
            local_scheduler = False
        else:
            local_scheduler = True

        sl.run(
            local_scheduler=local_scheduler,
            main_task_cls=sv_dada2.Workflow_DADA2,
            cmdline_args=[
                '--workers={}'.format(
                    args.workers
                    ),
                '--working-dir={}'.format(
                    args.working_dir
                    ),
                '--destination-dir={}'.format(
                    args.destination_dir
                    ),
                '--manifest={}'.format(
                    args.manifest
                    ),
            ]
        )

    def refpkg(self):
        """Make a reference package appropriate for pplacer or other pipelines."""
        args = self.__args__

        if args.luigi_manager:
            local_scheduler = False
        else:
            local_scheduler = True

        sl.run(
            local_scheduler=local_scheduler,
            main_task_cls=refpkg.WorkflowMakeRefpkg,
            cmdline_args=[
                '--sequence-variants-path={}'.format(
                    args.sequence_variants
                    ),
                '--repo-filtered-seq-info={}'.format(
                    args.repo_filtered_seq_info
                    ),
                '--repo-filtered={}'.format(
                    args.repo_filtered
                    ),
                '--repo-genomes-seq-info={}'.format(
                    args.repo_genomes_seq_info
                    ),
                '--repo-genomes={}'.format(
                    args.repo_genomes
                    ),
                '--new-refpkg-path={}'.format(
                    args.refpkg_destdir
                    ),
                '--new-refpkg-name={}'.format(
                    args.refpkg_name
                    ),
                '--working-dir={}'.format(
                    args.working_dir
                    ),
                '--min-id-genomes={}'.format(
                    args.min_id_genomes
                    ),
                '--min-id-filtered={}'.format(
                    args.min_id_filtered
                    ),
                '--min-id-unnamed={}'.format(
                    args.min_id_unnamed
                    ),
                '--min-best={}'.format(
                    args.min_best
                    ),
                '--workers={}'.format(
                    args.workers
                ),
            ]
        )

    def ncbi_16s(self):
        """Update a repository of 16S sequences from NCBI NT."""
        print("Starting NCBI_16s")
        args = self.__args__

        if args.luigi_manager:
            local_scheduler = False
        else:
            local_scheduler = True

        sl.run(
            local_scheduler=local_scheduler,
            main_task_cls=ncbi_16s.Workflow_NCBI_16s,
            cmdline_args=[
                '--ncbi-email={}'.format(
                    args.ncbi_email
                    ),
                '--repo-url={}'.format(
                    args.repo_secret
                ),
                '--example-seqs={}'.format(
                    args.example_seqs
                ),
                '--working-dir={}'.format(
                    args.working_dir
                    ),
                '--workers={}'.format(
                    args.workers
                ),
            ]
        )

    def placement(self):
        """Place sequence variants on a reference package."""
        args = self.__args__

        if args.luigi_manager:
            local_scheduler = False
        else:
            local_scheduler = True

        cmdline_args = [
                '--sv-fasta={}'.format(
                    args.sequence_variants
                    ),
                '--working-dir={}'.format(
                    args.working_dir
                    ),
                '--destination-dir={}'.format(
                    args.destination_dir
                    ),
                '--refpkg-tgz={}'.format(
                    args.refpkg_tgz
                    ),
                '--seq-map-csv={}'.format(
                    args.seq_map_csv,
                    ),
                '--workers={}'.format(
                    args.workers
                ),
            ]

        if args.sv_weights_csv:
            cmdline_args.append(
                '--sv-weights-csv={}'.format(
                    args.sv_weights_csv,
                    ),
            )

        sl.run(
            local_scheduler=local_scheduler,
            main_task_cls=placement.Workflow_Placement,
            cmdline_args=cmdline_args
        )

    def classify(self):
        """Classify sequence variants using a placement."""
        args = self.__args__

        if args.luigi_manager:
            local_scheduler = False
        else:
            local_scheduler = True

        cmdline_args = [
                '--sv-fasta={}'.format(
                    args.sequence_variants
                    ),
                '--working-dir={}'.format(
                    args.working_dir
                    ),
                '--jplace={}'.format(
                    args.jplace
                ),
                '--destination-dir={}'.format(
                    args.destination_dir
                    ),
                '--refpkg-tgz={}'.format(
                    args.refpkg_tgz
                    ),
                '--seq-map-csv={}'.format(
                    args.seq_map_csv,
                    ),
                '--workers={}'.format(
                    args.workers
                ),
            ]

        if args.sv_weights_csv:
            cmdline_args.append(
                '--sv-weights-csv={}'.format(
                    args.sv_weights_csv,
                    ),
            )

        if args.labels:
            cmdline_args.append(
                '--labels={}'.format(
                    args.labels,
                    ),
            )

        sl.run(
            local_scheduler=local_scheduler,
            main_task_cls=classify.Workflow_Classify,
            cmdline_args=cmdline_args
        )

def main():
    """Entrypoint for main script."""
    MALIAMPI()


if __name__ == "__main__":
    main()
