#!/usr/bin/env python3

import re
import os
import sys
import gzip
import argparse
from argparse import RawTextHelpFormatter


def is_gz_file(filepath):
    """Determine whether the provided file is GZipped or not.
    Returns 'rt' if GZipped file
    else returns 'r' """
    with gzip.open(filepath, 'r') as fh:
        try:
            # Gzip file
            fh.read(1)
            return 'rt'
        except gzip.BadGzipFile:
            # regular file
            return 'r'


def parse_pfam_res(infile):
    """Get all Pfams for each protein"""
    # Handle regular and gzip file
    read_mode = is_gz_file(infile)
    if read_mode == 'rt':
        fi = gzip.open(infile, read_mode)
    else:
        fi = open(infile, read_mode, encoding = 'utf-8')

    res = {}  # store the result => {prot: [pfams]}

    # Browse the file
    for line in fi.readlines():
        if line.startswith('#'):
            continue

        li = line.rstrip().split(None)  # None allows to discard white spaces!!

        prot_id = li[0]
        pfam_acc = li[3].split('.')[0]  # keep PF00389, not PF00389.33

        # Save the result
        try:
            res[prot_id] += [pfam_acc]
        except KeyError:
            res[prot_id] = [pfam_acc]

    fi.close()

    return res


def print_results(res, outfile):
    """A function to format the results stored as dict"""
    with open(outfile, 'w', encoding = 'utf-8') as fo:
        fo.write('\t'.join(['protein', 'Pfams']) + '\n')

        for prot, pfams in res.items():
            fo.write(prot + '\t' + ','.join(pfams) + '\n')

    return True


if __name__ == '__main__':
    __description__ = 'Read an annotation file produced by HMMsearch and ' \
                      'outputs a TSV with 2 columns:\n' \
                      'sequence ID <tab> Pfam domain(s)\n'

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description = __description__)
    parser.add_argument('-i', '--input', help="Hmmsearch results",
                        required = True, metavar = '')
    parser.add_argument('-o', '--output', required = True, metavar = "")

    args = parser.parse_args()

    try:
        # Checks
        if os.path.exists(args.output):
            raise FileExistsError('\nThe name provided for the output file ' +
                                'exists. Please remove existing file or ' +
                                 'provide another name\n')
        # Parse
        results = parse_pfam_res(args.input)

        # Print
        print_results(results, args.output)

    except Exception as e:
        print(e)
        sys.exit(1)

