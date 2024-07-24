# =========================================== Modules ============================================ #

import os
import time
from Bio import SeqIO


# ================================================================================================ #

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


def detect_files(files):
    """
    Checks if files is a file containing data or just containing path to other files
    """
    fs = None
    for file in files:
        if determine_file(file):
            fs = get_files_from_arg(file)

    if fs:
        return fs
    else:
        return files


def read_fasta_files(files, output):
    """
    Reads a fasta files and concatenate their content in a unique file
    """
    i = 1
    with open(output[0], "w") as out:
        with open(output[1], "w") as index:
            for file in files:
                records = SeqIO.parse(open(file), "fasta")
                name = "-".join(file.split("-")[1:2]).split("/")[-1]
                name = name.replace("_", "-") + "-1"
                for record in records:
                    if "METDB" in file:
                        if "TRINITY" in record.id:
                            record.id = record.id.replace("TRINITY", name.upper())
                            record.description = record.description.replace("TRINITY", name.upper())
                        record.id = record.id.replace("_", "-")
                        record.description = record.description.replace("_", "-")
                    if len(record.id.split("|")) > 1:
                        record.id = record.id.split("|")[0]
                    index.write(f"{i}\t{record.id}\t{record.description}\n")
                    record.id = str(i)
                    record.description = ""
                    SeqIO.write(record, out, "fasta")
                    i += 1


def main():
    """
    Main program function
    """
    with open(str(snakemake.log), "w") as log:
        s = time.time()
        log.write("*** Getting input and output files ***\n")
        files = detect_files(snakemake.input)
        output = list(snakemake.output)
        log.write("*** Concatenation ***\n")
        read_fasta_files(files, output)
        e = time.time()
        log.write(f"Concatenation of your files done in {round(e - s, 2)} seconds")


if __name__ == '__main__':
    main()
