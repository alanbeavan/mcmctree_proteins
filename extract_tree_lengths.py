#!/usr/bin/env python3
"""Docstring."""

import ete3
import glob

def main():
    """Do the things."""
    table = ["genefam,dist"]
    treefiles = glob.glob("Pinniped_phylogeny_files/gene_trees/trees/*.treefile")
    for tf in treefiles:
        og = tf.split("/")[3].split(".")[0]
        tree = ete3.Tree(tf)
        tree.set_outgroup(tree.get_midpoint_outgroup())
        farthest, dist = tree.get_farthest_node()
        table.append(og + "," + str(dist*2))
    with open("tree_lengths.csv", "w", encoding = "utf8") as out:
        out.write("\n".join(table) + "\n")
        

if __name__ == "__main__":
    main()
