# =========================================== Modules ============================================ #

import csv
import io
import gzip
import time


# ================================================================================================ #

def percentage(part, whole):
    """
    Retrieves the percentage
    """
    percent = 100 * float(part) / float(whole)
    return str(round(percent, 2)) + "%"


def get_eval(row):
    """
    Retrieves the E-value exponent
    """
    if type(row) == dict:
        e = row["evalue"]
    else:
        e = row
    if e == "0.0":
        return 0
    return int(e.split("e")[1])


def filter(inputfile, outputfile, cov, ident, eval):
    """
    Filtration Function
    """
    al_ssn, al_filt, nb_nssn, nb_nfilt = 0, 0, 0, 0
    n_ssn = set()
    n_filt = set()
    snu = set()
    ev = int(eval.split("e")[1])

    fieldnames = ["qseqid", "qlen", "qstart", "qend", "sseqid", "slen",
                  "sstart", "send", "length", "pident", "ppos", "score", "evalue",
                  "bitscore"]

    if inputfile.endswith(".gz"):
        i = gzip.open(inputfile, mode="rt")
        reader = csv.DictReader(i, delimiter="\t",
                                fieldnames=fieldnames)
    else:
        reader = csv.DictReader(open(inputfile), delimiter="\t",
                                fieldnames=fieldnames)

    writer = csv.DictWriter(open(outputfile, "w"), delimiter="\t",
                            fieldnames=fieldnames)

    for row in reader:

        al_ssn += 1
        n = [row["qseqid"], row["sseqid"]]
        e_row = get_eval(row)
        ns = sorted(n)
        nu = ns[0]+"="+ns[1]

        if n[0] not in n_ssn:
            n_ssn.add(n[0])
            nb_nssn += 1

        if (
                row["pident"]
                and row["ppos"]
                and e_row <= ev
                and float(row["pident"]) >= ident
                and float(row["ppos"]) >= cov
                and n[0] != n[1]
                and nu not in snu
        ):
            writer.writerow(row)
            al_filt += 1
            snu.add(nu)
            if n[0] not in n_filt:
                n_filt.add(n[0])
                nb_nfilt += 1

    return al_ssn, al_filt, nb_nssn, nb_nfilt


def filter_file(identity, overlap, eval, ssn, log, filtr, stats):
    """
    Filters the File
    """

    for i in identity:
        for j in overlap:
            for h in eval:
                print("*** FILTERING FILE ***", file=log)
                print(f"filtering infos ||| coverage : {j}%, identity : {i}%, eval : {h}", file=log)
                output = filtr[0]
                filtr.remove(output)
                al_ssn, al_filt, nb_nssn, nb_nfilt = filter(ssn, output, j, i, h)
                rm = al_ssn - al_filt
                rmnb = nb_nssn - nb_nfilt
                print(f"nb of alignments in base SSN : {al_ssn}", file=log)
                print(f"nb of alignments in filtered SSN : "
                      f"{al_filt} ({percentage(al_filt, al_ssn)})",
                      file=log)
                print(f"nb of alignments removed : {rm} ({percentage(rm, al_ssn)})", file=log)
                print(f"nb of nodes in base SSN : {nb_nssn}", file=log)
                print(f"nb of nodes in filtered SSN : {nb_nfilt} ({percentage(nb_nfilt, nb_nssn)})",
                      file=log)
                print(f"nb of nodes removed : {rmnb} ({percentage(rmnb, nb_nssn)})", file=log)
                print("*** Saving Stats ***\n", file=log)
                output_stats = stats[0]
                stats.remove(output_stats)
                save_stats(output_stats, al_ssn, al_filt, nb_nssn, nb_nfilt)


def save_stats(output, al_ssn, al_filt, nb_nssn, nb_nfilt):
    """
    Saves filtration stats
    """
    rm = al_ssn - al_filt
    rmnb = nb_nssn - nb_nfilt
    with open(output, "w") as f:
        f.write(f"nb of alignments in base SSN : {al_ssn}\n")
        f.write(f"nb of alignments in filtered SSN : {al_filt} ({percentage(al_filt, al_ssn)})\n")
        f.write(f"nb of alignments removed : {rm} ({percentage(rm, al_ssn)})\n")
        f.write(f"nb of nodes in base SSN : {nb_nssn}\n")
        f.write(f"nb of nodes in filtered SSN : {nb_nfilt} ({percentage(nb_nfilt, nb_nssn)})\n")
        f.write(f"nb of nodes removed : {rmnb} ({percentage(rmnb, nb_nssn)})\n")


# def diamond2graph(diamond_output):
#    ids = set()
#    fieldnames = "from;to;query_length;subject_length;alignment_len;pident;evalue;bitscore"
#    edges = open(f"{diamond_output}.edges", "w")
#    vertices = open(f"{diamond_output}.vertices", "w")
#    edges.write(f"{fieldnames}\n")
#    vertices.write("name;prefix\n")
#    for line in read_file(diamond_output):
#        tab = line.split("\t")
#        qsid, qlen, ssid, slen, leng, pid, eval, bscore = tab[0], tab[1], tab[4],\
#                                                          tab[5], tab[8], tab[9],\
#                                                          tab[12], tab[13]
#        edges.write(f"{qsid};{ssid};{qlen};{slen};{leng};{pid};{eval};{bscore}")

#        if qsid not in ids:
#            ids.add(qsid)
#        if ssid not in ids:
#            ids.add(ssid)

#   ids = sorted(ids)

#    for i in ids:
#        prefix = i.split(".")[0]
#        vertices.write(f"{i};{prefix}\n")


def main():
    """
    Main program function
    """
    with open(str(snakemake.log), "w") as log:
        s = time.time()
        log.write("*** Filtering information provided ***\n")
        filter_file(snakemake.params.identity,
                    snakemake.params.overlap,
                    snakemake.params.evalue,
                    str(snakemake.input), log,
                    snakemake.output.filtr,
                    snakemake.output.stats)
        e = time.time()
        log.write(f"Operations done in {round(e - s, 2)} seconds")


if __name__ == '__main__':
    main()
