#!/usr/bin/env python
import argparse
import logging
import sys
import luigi
import sciluigi as sl

from subcommands import refpkg, ncbi_16s

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

            # refpkg command
            subparser_refpkg = self.__subparsers__.add_parser(
                'refpkg',
                description="""Make a reference package
                                appropriate for pplacer or other pipelines.""")
            refpkg.build_args(subparser_refpkg)

            # repo_16s command
            subparser_ncbi_16s = self.__subparsers__.add_parser(
                'ncbi_16s',
                description="""Update a repository of 16S sequences from NCBI NT.
                """)
            ncbi_16s.build_args(subparser_ncbi_16s)

            self.__args__ = parser.parse_args()
            if self.__args__.command == 'refpkg':
                self.refpkg()
            if self.__args__.command == 'ncbi_16s':
                self.ncbi_16s()

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
                '--repo-seq-info={}'.format(
                    args.repo_seq_info
                    ),
                '--repo-seqs-filtered={}'.format(
                    args.repo_seqs_filtered
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
                '--min-id-types={}'.format(
                    args.min_id_types
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
                '--working-dir={}'.format(
                    args.working_dir
                    ),
                '--workers={}'.format(
                    args.workers
                ),
            ]
        )


def main():
    """Entrypoint for main script."""
    MALIAMPI()


if __name__ == "__main__":
    main()
