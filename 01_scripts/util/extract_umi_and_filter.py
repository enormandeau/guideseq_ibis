#!/usr/bin/env python3
"""Extract UMI+ info from PE reads and filter

Usage:
    <program> fqz1 fqz2 trim_length output_fasta
"""

# Modules
import gzip
import sys

# Classes
class Fasta(object):
    """Fasta object with name and sequence
    """

    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def write_to_file(self, handle):
        handle.write(">" + self.name + "\n")
        handle.write(self.sequence + "\n")

    def __repr__(self):
        return self.name + " " + self.sequence[:31]

class Fastq(object):
    """Fastq object with name, sequence, name2, and quality string
    """

    def __init__(self, name, sequence, name2, quality):
        self.name = name
        self.sequence = sequence
        self.name2 = name2
        self.quality = quality

    def getShortname(self, separator):
        if separator:
            self.temp = self.name.split(separator)
            del(self.temp[-1])
            return separator.join(self.temp)

        else:
            return self.name

    def write_to_file(self, handle):
        handle.write(self.name + "\n")
        handle.write(self.sequence + "\n")
        handle.write(self.name2 + "\n")
        handle.write(self.quality + "\n")

    def __repr__(self):
        return self.name + " " + self.sequence[:31]

# Defining functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)

def fastq_iterator(infile):
    """Takes a fastq file infile and returns a fastq object iterator

    Requires fastq file with four lines per sequence and no blank lines.
    """
    
    with myopen(infile) as f:
        while True:
            name = f.readline().strip()

            if not name:
                break

            seq = f.readline().strip()
            name2 = f.readline().strip()
            qual = f.readline().strip()
            yield Fastq(name, seq, name2, qual)

# Expected machinery
## Read 1
umi = ""
alien = "CCATCTCATCCCTGC"

## Read 2
odn_plus = "CGTTATTAACATATGACAACTCAATTAAAC"
odn_minus = "TTGAGTTGTCATATGTTAATAACGGTAT"

# Misc
n = 0
odn_len = 0
num_kept = 0
num_flushed = 0

# Parse user input
try:
    fqz1 = sys.argv[1]
    fqz2 = sys.argv[2]
    trim_length = int(sys.argv[3])
    output_fasta = sys.argv[4]
except:
    print(__doc__)
    sys.exit(1)

fq1 = fastq_iterator(fqz1)
fq2 = fastq_iterator(fqz2)
outfile = myopen(output_fasta, "wt")

for pair in zip(fq1, fq2):
    s1, s2 = pair
    keep = True
    names = [x.name.split(" ")[0] for x in pair]

    assert names[0] == names[1], f"\n\nGUIDEseq IBIS: Sequences in fastq files are not synced\
\n\n{fqz1}: {names[0]}\n{fqz2}: {names[1]}"

    # TODO filter if ALIEN or dsODN found more than once
    # TODO permit diffs in ALIEN and dsODN

    if alien in s1.sequence[8:24]:
        #print("Alien found")
        if odn_plus in s2.sequence[:len(odn_plus) + 1]:
            #print("odn+ found")
            num_kept += 1
            odn_len = len(odn_plus)

        elif odn_minus in s2.sequence[:len(odn_minus) + 1]:
            #print("odn- found")
            num_kept += 1
            odn_len = len(odn_minus)

        else:
            #print("odn not found")
            num_flushed += 1
            keep = False

    else:
        #print("Alien not found")
        num_flushed += 1
        keep = False

    if keep:
        umi = s1.sequence[:8]
        first8 = s1.sequence[24:32]
        sequence = s2.sequence[odn_len+1: ]

        if len(sequence) < trim_length:
            #print("seq too short")
            num_kept -= 1
            continue

        # Write fasta file
        n += 1
        s = Fasta(f"seq{n}_{umi}_{first8}", sequence[:trim_length])
        s.write_to_file(outfile)

percentage = round((100 * num_kept) / (num_kept + num_flushed), 2)
print(f"Retained {percentage}% ({num_kept} / {num_kept + num_flushed})")
