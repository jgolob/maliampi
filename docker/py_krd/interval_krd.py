#!/usr/bin/env python3
import argparse
import os
import subprocess
import csv

# ------------------------------------------------------------------
#      Take a set of NEW specimen-jplace files and relate to extant
#      specimen-jplace files. N:1 (i.e. for each new relate to each old
#      but not old to old or new to new.)
# ------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description='Integrate a new jplace into extant jplace')
    parser.add_argument(
        '--new-jplace',
        required=True
    )
    parser.add_argument(
        '--old-jplaces',
        required=True
    )

    args = parser.parse_args()

    new_jplace = args.new_jplace

    #
    old_jplaces = [
        jp for jp in args.old_jplaces.split()
        if jp.endswith('jplace.gz')
    ]

    successful_pairs = []
    for ojp_i, ojp_fn in enumerate(old_jplaces):
        sp = subprocess.run([
            'gappa',
            'analyze',
            'krd',
            '--file-prefix', "v{:05d}__".format(ojp_i),
            '--matrix-format', 'list',
            '--jplace-path', new_jplace,
        ] + old_jplaces)
        if (sp.returncode == 0):
            successful_pairs.append(
                'v{:05d}__krd_matrix.csv'.format(ojp_i)
            )
    # Great! We are now completed. Extract out the successful pairs
    distance_long = []
    for d_fn in successful_pairs:
        with open(d_fn, 'rt') as in_h:
            for r in csv.reader(in_h, delimiter='\t'):
                if r[0] != r[1]:
                    distance_long.append({
                        'specimen_1': r[0],
                        'specimen_2': r[1],
                        'krd': float(r[2])
                    })
    # Finally output
    with open("v{}.krd_long.csv".format(os.path.basename(new_jplace)).replace('.jplace.gz', ''), 'wt') as out_h:
        w = csv.DictWriter(out_h, ['specimen_1', 'specimen_2', 'krd'])
        w.writeheader()
        w.writerows(distance_long)


if __name__ == "__main__":
    main()
