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
    fragment_length = int(sys.argv[2])
    output_file = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Find DSBs
min_cov = 10
window_size = 10

# Load gene annotation file
genes = [x.strip().split("\t") for x in
        open("00_archive/human_GRCh38.p14_genes_simplified.tsv").readlines()]

genes_dict = defaultdict(lambda: defaultdict(list))
bin_size = 10000

for g in genes:
    if g[0] == ("#Name"):
        continue

    chrom, _from, _to, gene = g[1], int(g[2]), int(g[3]), g[6]

    bins = range(_from // bin_size, (_to // bin_size)+1)

    for b in bins:
        genes_dict[chrom][b].append((_from, _to, gene))

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
    outfile.write("0Sample\tChromosome\tSite1\tSite2\tDist\tCnt1\tCnt2\tRatio\tGene\n")

    for s in sorted(sites):
        chrom, pos = s
        count = sites[s]
        sites.pop(s)

        # Accept only positions close to the expected distance. The fragments
        # should align one fragment-length apart for it to be a true DSB. We
        # can accept an error of 10-20bp
        error = 10

        for p in range(pos+fragment_length-error, pos+fragment_length+error+1):
            _id = (s[0], p)

            if _id in sites:
                num_sites += 1

                # Get closest gene
                gene_set = set()
                break_pos = _id[1]

                # Add overlapping gene info
                for gene in genes_dict[chrom][break_pos // bin_size]:
                    if break_pos in range(gene[0], gene[1]+1):
                        gene_set.add(gene[-1])

                gene_name = ",".join(sorted(list(gene_set)))
                if not gene_name:
                    gene_name = "-na-"

                read_ratio = sites[_id] / count
                read_ratio = round(read_ratio, 1)

                # Write to file
                outfile.write("\t".join([str(x) for x in [sample, s[0], s[1], _id[1],
                    _id[1] - s[1] - fragment_length, count, sites[_id], read_ratio, gene_name]]) + "\n")

    # Report stats for file
    print(f"Found {num_sites} sites for {input_file.split('/')[1].split('_')[0]}")
