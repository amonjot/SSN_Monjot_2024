# =========================================== Modules ============================================ #

import os
import csv
import time


# ================================================================================================ #

def read_csv(file):
    """
    Reads a file as CSV
    """
    with open(file, "r") as f_in:
        reader = csv.DictReader(f_in, delimiter="\t")
        yield from reader


def adapt_row(row, columns, i_dict):
    """
    Adapt incomplete row with columns header
    """

    if row[columns[0]] not in i_dict:
        return False, False, False
    n = i_dict[row[columns[0]]]
    c = columns[1:]
    d = {k: v for k, v in row.items() if k in c}

    for k, v in d.items():
        if v is None:
            d[k] = "NA"

    return d, n, c


def determine_file(file):
    """
    Checks if the file contains data or a list of path to files
    """
    if type(file) != str:
        file = str(file)
    with open(file) as f:
        line = f.readline()
        return bool(os.path.exists(line.strip()))


def get_files_from_arg(file):
    """
    Retrieves a list of files
    """
    if type(file) != str:
        file = str(file)
    files = []
    with open(file, "r") as f_in:
        for line in f_in:
            files.append(line.strip())

    return files


def treat_annotation(an_files, columns, indices):
    """
    Treats annotation files
    """
    at_dict = {}
    i_dict = {}

    with open(indices, "r") as f:
        for line in f:
            llist = line.split("\t")[:2]
            i_dict[llist[1]] = llist[0]

    for file in an_files:
        n = file.split("/")[-1].split(".")[0]
        for row in read_csv(file):
            d, nc, c = adapt_row(row, columns, i_dict)
            if d:
                if n not in at_dict.keys():
                    at_dict[n] = {}
                if nc not in at_dict[n].keys():
                    at_dict[n][nc] = {k: set() for k in d.keys()}
                for k, v in d.items():
                    at_dict[n][nc][k].add(v)

    return at_dict


def create_attributes_dict(an_files, columns, indices):
    """
    Creates a dictionary containing all attributes
    """
    if len(an_files) != 1 or not determine_file(an_files):
        return treat_annotation(an_files, columns, indices)
    files = get_files_from_arg(an_files)
    return treat_annotation(files, columns, indices)


def search_output(k, output):
    """
    Retrieves output corresponding to k
    """
    for o in output:
        if k in o:
            return o


def fill_tmp_dict(k, v, columns):
    d = {}
    for key, val in v.items():
        val = "".join(val) if len(val) == 1 else ",".join(val)
        d[key] = val

    for c in columns:
        if c not in d.keys() and k not in d.values():
            d[c] = k
    return d


def save_attributes(at_dict, columns, outputs):
    """
    Saves attributes in a file
    """
    path = "./results/attributes/"
    if not os.path.exists(path):
        os.mkdir(path)

    for k, v in at_dict.items():
        output = search_output(k, outputs)
        with open(output, "w") as f:
            writer = csv.DictWriter(f, delimiter=";", fieldnames=columns)
            writer.writeheader()
            for k2, v2 in v.items():
                tmp = fill_tmp_dict(k2, v2, columns)
                writer.writerow(tmp)


def main():
    """
    Main program function
    """
    with open(str(snakemake.log), "w") as log:
        s = time.time()
        log.write("*** Getting input files and parameters ***\n")
        columns = snakemake.params.columns
        files = snakemake.input.an_files
        indices = snakemake.input.indices
        log.write("*** Getting Attributes from files ***\n")
        at_dict = create_attributes_dict(files, columns, indices)
        log.write("*** Saving Attributes ***\n")
        save_attributes(at_dict, columns, snakemake.output)
        e = time.time()
        log.write(f"Operations done in {round(e - s, 2)} seconds")


if __name__ == '__main__':
    main()
