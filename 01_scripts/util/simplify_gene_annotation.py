#!/usr/bin/env python3
"""Simplify gene annotation from human_GRCh38.p14_genes_simplified.tsv

Usage:
    <program> input_file output_file

Where input_file is human_GRCh38.p14_genes_simplified.tsv
"""

# Module
import sys

# Functions
def overlaps(r1, r2):
    """Return if range 1 (r1 = [from<int>, to<int>]) overlaps range 2 (r2 = [from<int>, to<int>])
    """
    return (r1[1] >= r2[0]) & (r2[1] >= r1[0])

# Parse user input
try:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Simplify
prev = ("", "" ,"", "", "", "", "")

with open(output_file, "wt") as outfile:
    with open(input_file) as infile:
        for line in infile:

            l = line.strip().split("\t")

            # Header
            if line.startswith("#"):
                outfile.write(line)
                continue

            # Confirm chrom, pos, gene are the same
            if (l[0], l[1], l[6]) == (prev[0], prev[1], prev[6]):

                # Confirm the regions overlap and grow the region if necessary
                if overlaps((l[2], l[3]), (prev[2], prev[3])):
                    prev[2] = min([l[2], prev[2]])
                    prev[3] = max([l[3], prev[3]])

            # Write to file
            else:
                if prev[0]:
                    outfile.write("\t".join(prev) + "\n")

                prev = l

    outfile.write("\t".join(prev) + "\n")
