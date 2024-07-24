def read_file(file):
    """
    Reads a file and retrieves its content line by line
    """
    with open(file, "r") as f:
        for line in f:
            yield line.strip()


def diamond2graph(diamond_output):
    """
    Converts Diamond output into igraph-python input
    """
    fieldnames = "from;to;query_length;subject_length;alignment_len;pident;evalue;bitscore"
    for file in diamond_output:
        ids = set()
        edges = open(f"{file}.edges", "w")
        vertices = open(f"{file}.vertices", "w")
        edges.write(f"{fieldnames}\n")
        vertices.write("name\n")

        for line in read_file(file):
            tab = line.split("\t")
            qsid, qlen, ssid, slen, leng, pid, eva, bscore = tab[0], tab[1], tab[4], \
                                                             tab[5], tab[8], tab[9], \
                                                             tab[12], tab[13]
            edges.write(f"{qsid};{ssid};{qlen};{slen};{leng};{pid};{eva};{bscore}\n")

            if qsid not in ids:
                ids.add(qsid)
            if ssid not in ids:
                ids.add(ssid)

        ids = sorted(ids)

        for i in ids:
            vertices.write(f"{i}\n")


def main():
    """
    Main program function
    """
    diamond2graph(snakemake.input)


if __name__ == '__main__':
    main()
