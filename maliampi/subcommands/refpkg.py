#!/usr/bin/python


def build_args(parser):
    parser.add_argument(
        '--sequence-variants',
        help="""Path to sequence variants (in FASTA format)
            for which we need a reference set created""",
        required=True
        )
    parser.add_argument(
        '--repo-seqs-filtered',
        help="""Path(s) to repository sequences with trusted annotations
            from which we should recruit. FASTA format expected""",
        type=str,
        nargs='+',
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
        '--working_dir',
        help="""Path of a suitable working directory
        (defaults to the current working directory)"""
    )
    parser.add_argument(
        '--min_id_types',
        default=0.8,
        type=float,
        help='Min percent identity when recruiting from type strain 16S'
    )
    parser.add_argument(
        '--min_id_filtered',
        default=0.8,
        type=float,
        help='Min percent identity when recruiting from filtered 16S'
    )
    parser.add_argument(
        '--min_id_unnamed',
        default=0.97,
        type=float,
        help='Min percent identity when recruiting from unnamed / no taxonomy 16S'
    )
    parser.add_argument(
        '--min_best',
        default=1.0,
        type=float,
        help='Min percent identity to consider a query matched'
    )

    #
    # parser.add_argument(
    #   '--taxdump',
    #   help='Path to taxdump.tar.gz to use',
    #   default = None
    # )
    #


def main():
    """Entrypoint for main script."""
    print('Hello!')


if __name__ == "__main__":
    main()
