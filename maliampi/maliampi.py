#!/usr/bin/python
import argparse
import logging
import sys

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
        print(self.__args__)


def main():
    """Entrypoint for main script."""
    MALIAMPI()


if __name__ == "__main__":
    main()
