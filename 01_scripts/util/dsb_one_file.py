#!/usr/bin/env python3
"""Find Double Stranded Breaks (DSBs)

Usage:
    <program> input_file output_file
"""

# Modules
# import gzip
from collections import defaultdict
import sys

# Parse input
try:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Find DSBs
min_cov = 20
window_size = 10

#TODO Load value from config_file
neighbour_size = 200
coverages = defaultdict(int)
data = [x.strip().split("\t") for x in open(input_file).readlines()]

for d in data:
    _id = (d[0], int(d[1]))
    coverages[_id] += 1

sorted_coverages = sorted([(x[1], x[0]) for x in coverages.items()], reverse=True)

# Get site with higest coverage
# Aggregate all the counts withing `window_size` bp on each side
# Remove site and neighbours from `coverages`
# Count this as a potential site and go down the list until you hit
# a site with less than `min_cov` reads
sites = {}
visited = set()

for site in sorted_coverages:
    cov = site[0]
    if cov < min_cov:
        break

    chrom, pos = site[1]

    if (chrom, pos) in visited:
        continue

    sites[(chrom, pos)] = 0

    for p in range(pos-window_size, pos+window_size+1):
        visited.add((chrom, p))

        if (chrom, p) in coverages:
            sites[(chrom, pos)] += coverages.pop((chrom, p), 0)

# Get DSB duets: pairs of sites with good coverages at ~100 bp one from each other
duets = 0
sample = input_file.split("/")[1].split("_")[0]
num_sites = 0

with open(output_file, "wt") as outfile:
    outfile.write("Sample\tChrom\tSite1\tSite2\tDistance\tCount1\tCount2\n")

    for s in sorted(sites):
        chrom, pos = s
        count = sites[s]
        sites.pop(s)

        # TODO Accept only positions close to the wanted distance
        # as opposted to smaller than or equal to a max distance
        for p in range(pos-neighbour_size, pos+neighbour_size+1):
            _id = (s[0], p)

            if _id in sites:
                num_sites += 1
                outfile.write("\t".join([str(x) for x in [sample, s[0], s[1], _id[1], _id[1] - s[1], count, sites[_id]]]) + "\n")

    print(f"Found {num_sites} sites for {input_file.split('/')[1].split('_')[0]}")
