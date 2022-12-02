#!/bin/bash
# Map fasta reads to genome

# Global variables
BASE=$(basename $1)
NCPUS="$2"
GENOME_FILE="$3"
GENOME_FOLDER="03_genome"
EXTRACTED_FOLDER="06_extracted"
ALIGNED_FOLDER="07_aligned"

# BWA
bwa mem -t "$NCPUS" -T 10 \
    "$GENOME_FOLDER"/"$GENOME_FILE" \
    "$EXTRACTED_FOLDER"/"$BASE".fasta.gz |
    samtools view -S -q 1 -F 4 -F 256 -F 2048 - |
    sort -k 3,4 -V |
    grep -vP "^@" \
    > "$ALIGNED_FOLDER"/"${BASE%.fasta.gz}".sam
