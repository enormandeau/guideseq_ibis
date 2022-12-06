#!/usr/bin/env python3
"""Remove GUIDE-seq PCR duplicates from SAM alignment

Usage:
    <program> input_file output_file
"""

# Modules
# import gzip
import sys

# Parse input
try:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Dedup
data = [x.strip().split("\t") for x in open(input_file)]
data = sorted([(x[2], x[3], "_".join(x[0].split("_")[1:]), x[9]) for x in data])

# TODO Permit errors
deduplicated = sorted(list(set(data)))

print(f"{input_file} has {len(deduplicated)} ({round(100*len(deduplicated) / len(data), 2)}%) unique reads.")

with open(output_file, "wt") as outfile:
    for d in deduplicated:
        outfile.write("\t".join(d) + "\n")
