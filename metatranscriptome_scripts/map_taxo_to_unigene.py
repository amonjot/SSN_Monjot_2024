#!/usr/bin/env python3
# -*- encoding : utf-8 *-*

import os
import sys
import argparse
from argparse import RawTextHelpFormatter


RANK_d = { 'd': 'superkingdom', 'k': 'kingdom','p': 'phylum','c': 'class',
          'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species',
          '-': 'no rank'}


def parse_bed(infile):
    """A function to get the relation between protein and Unigene
    Returns a dict where key = unigene; value = list of proteins"""

    unigene2prot = {}

    with open(infile, 'r', encoding='utf-8') as fi:
        for line in fi.readlines():
            li = line.rstrip().split("\t")

            if len(li) <= 2:
                # Header line or something else
                continue

            unigene = li[0]
            protein = li[3].split(';')[0].replace("ID=", "")

            try:
                unigene2prot[unigene] += [protein]
            except KeyError:
                unigene2prot[unigene] = [protein]

    return unigene2prot


def parse_taxonomy_table(infile):
    """A function that parses the taxonomy table provided by MetaEuk
    Input format = TSV
    protein ID; taxid; best rank; best rank name; NCBI Taxonomy

    prot_000098547 98059 genus Dinobryon -_cellular organisms;d_Eukaryota;
                                         -_Stramenopiles;-_Ochrophyta;
                                         c_Chrysophyceae;o_Chromulinales;
                                         f_Dinobryaceae;g_Dinobryon

    Returns a dict with protein as key, value is a list"""

    prot_taxo = {}

    with open(infile, 'r', encoding='utf-8') as fi:
        for line in fi.readlines():
            li = line.rstrip().split("\t")
            if len(li) != 5:
                prot_taxo[li[0]] = [li[2], li[3], "NA"]
            else:
                prot_taxo[li[0]] = [li[2], li[3], li[4]]

    return prot_taxo


def associate_taxonomy_to_unigene(unigene2prot, prot_taxo):
    """This is the core function of the script. It takes a list of proteins
    per Unigene and a taxonomic affiliation of proteins. 
    It merge the affiliations to output a taxonomic affiliation PER Unigene"""

    # Declare the structure to store the results. A list of lists of size 2
    # with [Unigene, taxonomy]
    res = []

    for unigene, proteins in unigene2prot.items():
        if len(proteins) == 1:
            # One prot for the unigene, easy-peasy
            res.append([unigene,
                        '::'.join(prot_taxo[proteins[0]]) ])

        else:
            # Several proteins, aouch
            # 1. recover all taxo except 'no rank'
            taxos = []

            for protein in proteins:
                # Concatenate all fields associated to the prot' taxo
                concat_fields = '::'.join(prot_taxo[protein])

                # Gather all taxonomies together, work with the full taxo
                if concat_fields != "no rank::unclassified::NA":
                    taxos.append(prot_taxo[protein][2])

            # Case all Unclassified
            if len(taxos) == 0:
                res.append([unigene, concat_fields])

            # Case one unclassified and one with a taxo
            elif len(taxos) == 1:
                # Unclassified(s) + 1 classified
                res.append([unigene, concat_fields])

            # Case several proteins with an associated taxonomy
            else:
                # print(unigene)
                # 1. find the taxo with the shortest number of field
                shortest = ''
                size = 100000000
                for tax in taxos:
                    if len(tax.split(';')) < size:
                        shortest = tax

                # 2. browse the ranks of the worst assignation
                # Use the ranks present in the shortest taxonomy and test their
                # presence in the other taxonomy presents for this Unigene.
                # It is a way to do a 'LCA' using only the text
                lca_rank = ''  # variable to store the LCA

                for rank in shortest.split('::')[-1].split(";"):
                    # Loop over the taxonomies available

                    record = True  # Bool to record the current rank
                    for tax in taxos:
                        # Test presence of the current rank in the current taxo
                        if rank not in tax:
                            record = False

                    if record:
                        lca_rank += rank + ';'
                    else:
                        # useless to continue, so stop the loop
                        break

                lca_rank = lca_rank.rstrip(';')  # Clean the trailling ';'

                # 3. Get the value for the sub-fields 0 and 1, eg
                #  no rank::-_Stramenopiles or family::Thalassiosiraceae
                if ";" in lca_rank:
                    rank_id, rank_name = lca_rank.rsplit(";", maxsplit=1)[-1].split('_')
                else:
                    # for '-_cellular organisms' only
                    rank_id, rank_name = lca_rank.split('_')

                lca_rank = '::'.join([RANK_d[rank_id], rank_name, lca_rank])

                # 4. save
                res.append([unigene, lca_rank])
    return res


def print_taxonomy(result_list, outfile):
    """Simply print the results in a tabular file"""

    with open(outfile, 'w', encoding = 'utf-8') as fo:
        fo.write('\t'.join(['Unigene', 'Taxonomy']) + '\n')

        for result in result_list:
            fo.write('\t'.join(result) + '\n')

    return True


if __name__ == '__main__':
    __description__ = 'Infer a taxonomic affiliation to Unigenes from its ' \
                      'proteins\n'

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
                                     description = __description__)
    parser.add_argument('-i', '--input', help = "MetaEuk taxonomy table", \
                        required = True, metavar = '')
    parser.add_argument('-b', '--bed', help = "BED file provided by "\
                        "TransDecoder", required = True, metavar = '')
    parser.add_argument('-o', '--output', required = True, metavar = '')

    args = parser.parse_args()

    try:
        # Checks
        if os.path.exists(args.output):
            raise FileExistsError('\nThe name provided for the output file ' +
                                'exists. Please remove existing file or ' +
                                 'provide another name\n')

        # Get the relation between protein and Unigene
        unigene2prot = parse_bed(args.bed)

        # Get the taxonomic affiliation for the proteins
        prot_taxo = parse_taxonomy_table(args.input)

        # Merge all together
        res = associate_taxonomy_to_unigene(unigene2prot, prot_taxo)

        # print
        print_taxonomy(res, args.output)

    except Exception as e:
        print(e)
        sys.exit(1)
