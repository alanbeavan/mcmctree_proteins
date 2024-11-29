#!/usr/bin/env python3
"""Docstring."""

def get_file_data(filename):
    """stores the lines of the program with name filename as a list"""
    import sys
    try:
        with open(filename) as in_file:
            lines = []
            for line in in_file:
                lines.append(line.rstrip("\n"))
        return(lines)

    except IOError as e:
        print("{}\n Error opening {}. Terminating program.".format(e, filename))
        sys.exit(1)

def main():
    """Do the things."""
    stats = mod.get_file_data("clusters.csv")[1:]
    seqs = {}
    for line in stats:
        fields = line.split(",")
        ali = {}
        if fields[2] not in seqs:
            seqs[fields[2]] = read_fasta("Pinniped_phylogeny_files/alignments/" + fields[0] + ".faa")
        else:
            for key, value in read_fasta("Pinniped_phylogeny_files/alignments/" + fields[0] + ".faa").items():
                seqs[fields[2]][key] += value
    with open("partitioned_alignment.phy", "w", encoding = "utf8") as out:
        for key, value in seqs.items():
            out.write(str(len(value)) + "  " + str(len(value["cat"])) + "\n\n")
            for key1, value1 in value.items():
                #biggest key length is 22
                spaces = 24 - len(key1)
                out.write(key1)
                for i in range(spaces):
                    out.write(" ")
                out.write(value1 + "\n")
            out.write("\n\n")


if __name__ == "__main__":
    main()
