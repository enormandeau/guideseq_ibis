#!/usr/bin/env python3
"""Annotate GUIDE-seq hit sequences

Usage:
    <program> guideseq_ibis_report guide_infos genome output_file

guide_infos:
    - See example file in 02_infos
    - Potential header lines start with a pound (`#`) sign
    - has 3 tab-separated columns:
      - Sample
      - GuideName
      - GuideSequence
"""

# Modules
from collections import defaultdict
import gzip
import sys

# Classes
class Fasta(object):
    """Fasta object with name and sequence
    """

    def __init__(self, name, sequence):
        self.name = name
        self.shortname = self.name.split(" ")[0]
        self.sequence = sequence

    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")

    def __repr__(self):
        return self.shortname + " " + self.sequence[:31]

# Functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def fasta_iterator(input_file):
    """Takes a fasta file input_file and returns a fasta iterator
    """
    with myopen(input_file) as f:
        sequence = []
        name = ""
        begun = False

        for line in f:
            line = line.strip()

            if line.startswith(">"):
                if begun:
                    yield Fasta(name, "".join(sequence))

                name = line[1:]
                sequence = ""
                begun = True

            else:
                sequence += line

        if name != "":
            yield Fasta(name, "".join(sequence))

def hamming(s1, s2):
    """Find the Hamming distance btw. 2 strings. Substitutions only.
    
    TAKEN FROM: https://github.com/faircloth-lab/edittag
    ORIGINALLY FROM http://en.wikipedia.org/wiki/Hamming_distance
    
    """
    assert len(s1) == len(s2)
    return sum([ch1 != ch2 for ch1, ch2 in zip(s1, s2)])

def levenshtein(a, b):
    """Pure python version to compute the levenshtein distance between a and b.
    The Levenshtein distance includes insertions, deletions, substitutions; 
    unlike the Hamming distance, which is substitutions only.

    TAKEN FROM: https://github.com/faircloth-lab/edittag
    ORIGINALLY FROM:  http://hetland.org/coding/python/levenshtein.py

    """
    n, m = len(a), len(b)

    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a, b = b, a
        n, m = m, n

    current = range(n + 1)

    for i in range(1, m + 1):
        previous, current = current, [i] + [0] * n

        for j in range(1, n + 1):
            add, delete = previous[j] + 1, current[j - 1] + 1
            change = previous[j - 1]

            if a[j - 1] != b[i - 1]:
                change = change + 1

            current[j] = min(add, delete, change)

    return current[n]

def complement(s):
    new_s = []

    comp = dict(zip("ACGTN", "TGCAN"))

    for n in s.upper():
        new_s.append(comp[n])

    return "".join(new_s)

def revcomp(s):
    return complement(s)[::-1]

def find_best_position(target, query):
    best_pos = -1
    lowest_dist = len(query)

    for i in range(len(target) - len(query)):
        kmer = target[i: i+len(query)]
        dist = hamming(kmer, query)
        
        if dist < lowest_dist:
            best_pos = i
            lowest_dist = dist

    return (best_pos, lowest_dist)

# Parsing user input
try:
    guideseq_ibis_report = sys.argv[1]
    guide_infos = sys.argv[2]
    genome = sys.argv[3]
    output_file = sys.argv[4]

except:
    print(__doc__)
    sys.exit(1)

# Load hit infos
hits = [x.strip().split("\t") for x in open(guideseq_ibis_report).readlines()][1:]
hits_dict = defaultdict(set)

# Load guide infos
guides = [x.strip().split("\t") for x in open(guide_infos).readlines() if not x.startswith("#")]
guides = [x for x in guides if not x[2].lower() == "na"]
guides_dict = dict([(x[0], x[1: ]) for x in guides])

# Organize hits by chromosome to simplify sequence extraction and reduce memory
for h in hits:
    if h[0] in guides_dict:
        sample = h[0]
        chr_name = h[1]
        chr_id = h[2]
        pos = (int(h[3]), int(h[4]))

        hits_dict[chr_id].add((chr_name,) + pos + tuple(guides_dict[sample]))

# Document hits
annotations = dict()
sequences = fasta_iterator(genome)

for s in sequences:
    if s.shortname in hits_dict:

        for hit in hits_dict[s.shortname]:
            chr_id = s.shortname
            chr_name, _from, _to, guide_name, guide_seq = hit

            # Get best position in forward or reverse
            extract = s.sequence[_to-2*len(guide_seq)+2: _to+len(guide_seq)+2].upper()
            best_pos = find_best_position(extract, guide_seq)
            best_pos_rev = find_best_position(revcomp(extract), guide_seq)

            if best_pos[1] < best_pos_rev[1]:
                pos, dist = best_pos
                seq = extract[pos: pos+len(guide_seq)]
            else:
                pos, dist = best_pos_rev
                seq = revcomp(extract)[pos: pos+len(guide_seq)]

            # Create difference representation
            representation = []

            for i in range(len(seq)):
                if seq[i] == guide_seq[i]:
                    representation.append(seq[i])
                else:
                    representation.append("".join(["[", guide_seq[i], seq[i].lower(), "]"]))

            representation = "".join(representation)

            annotations[(chr_id, _from, _to)] = [guide_name, dist, guide_seq, representation]

            #print(" ".join([chr_name.rjust(13),
            #    str(_from).rjust(10),
            #    str(_to).rjust(10),
            #    str(dist).rjust(2),
            #    guide_seq,
            #    seq,
            #    representation]))

# Annotate report
lines = [x.strip().split("\t") for x in open(guideseq_ibis_report, "rt").readlines()]

with open(output_file, "wt") as outfile:
    for l in lines:
        if l[0] == "Sample":
            outfile.write("\t".join(l + ["GuideName", "Dist", "GuideSequence", "GuideAnnotation"]) + "\n")
            continue

        chr_id = l[2]
        _from = int(l[3])
        _to = int(l[4])

        key = (chr_id, _from, _to)

        if key in annotations:
            infos = [str(x) for x in annotations[key]]
        else:
            infos = ["na", "na", "na", "na"]

        outfile.write("\t".join(l + infos) + "\n")
