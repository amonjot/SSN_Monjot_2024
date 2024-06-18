#!/usr/bin/env python3

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


def parse_ko_res(infile):
    """Read and parse KoFamScan results
    
    As the are displayed by decreasing order of signicant results, this function
    checks:
    1. presence of '*'
    2. if not, take the first hit if evalue < 1e-5
    """
    # Read input
    read_mode = is_gz_file(infile)
    if read_mode == 'rt':
        fi = gzip.open(infile, read_mode)
    else:
        fi = open(infile, read_mode, encoding = 'utf-8')

    # Set-up some variables and structures
    res_list = []  # To store the results
    previous_prot_significative = ''  # Keep name of significative hit

    for line in fi.readlines():
        if line.startswith('#'):
            # Header line, so move to the next one
            continue

        li = line.rstrip().split("\t")

        if line[0] == "*":
            # Capture only the lines with positive match
            res_list.append([li[1], li[2], li[6].replace('"', ''),
                             "significant"])

            previous_prot_significative = li[1]

        elif li[1] == previous_prot_significative:
            # Next line as I do not want this protein
            continue
        else:
            # Take the first hit if evalue lower or equal to 1e-5
            if float(li[5]) <= 1e-5:
                res_list.append([li[1], li[2], li[6].replace('"', ''),
                                "best_hit"])

                previous_prot_significative = li[1]

    fi.close()
    return res_list


def print_results(res_list, outfile):
    """A function to format the results stored as list"""

    with open(outfile, 'w', encoding = 'utf-8') as fo:
        # print header
        fo.write('\t'.join(['protein', 'KO', 'KO_name', 'hit']) + '\n')

        for res in res_list:
            fo.write('\t'.join(res) + '\n')

    return True


if __name__ == '__main__':
    __description__ = 'Read an annotation file produced by KoFamScan and ' \
                      'outputs a TSV with 4 columns:\n' \
                      'sequence ID <tab> KO id <tab> KO name <tab> '\
                      'significative or best hit \n\n'\
                      '- Significative: score higher than the profile ' \
                      'threshold\n' \
                      '- Best hit: e-value <= 1e-5\n'

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description = __description__)
    parser.add_argument('-i', '--input', help = "KoFamScan result file, can " \
                        "be GZipped", required = True, metavar = '')
    parser.add_argument('-o', '--output', required = True, metavar = '')

    args = parser.parse_args()

    try:
        # Checks
        if os.path.exists(args.output):
            raise FileExistsError('\nThe name provided for the output file ' +
                                'exists. Please remove existing file or ' +
                                 'provide another name\n')

        # Parse
        results = parse_ko_res(args.input)

        # Print
        print_results(results, args.output)

    except Exception as e:
        print(e)
        sys.exit(1)
