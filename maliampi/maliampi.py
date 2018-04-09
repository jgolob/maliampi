#!/usr/bin/env python
import argparse
import logging
import sys
import luigi
import sciluigi as sl

from subcommands import refpkg


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

            self.__args__ = parser.parse_args()
            if self.__args__.command == 'refpkg':
                self.refpkg()

    def refpkg(self):
        """Make a reference package appropriate for pplacer or other pipelines."""
        args = self.__args__

        sl.run_local(
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


def main():
    """Entrypoint for main script."""
    MALIAMPI()


if __name__ == "__main__":
    main()
