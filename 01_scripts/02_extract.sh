#!/bin/bash
# Extract, filter, and annotate reads from ONE file

# Parameters
NCPUS=$1
MIN_LENGTH=$2

# Global variables
TRIMMED_FOLDER="05_trimmed"
EXTRACTED_FOLDER="06_extracted"

# Parallelize on all data files
ls -1 -S "$TRIMMED_FOLDER"/*_R1.fastq.gz |
    perl -pe 's/_R[12]\.fastq\.gz//' |
    parallel -k -j "$NCPUS" ./01_scripts/util/extract_umi_and_filter.py \
        {}_R{1,2}.fastq.gz "$MIN_LENGTH" "$EXTRACTED_FOLDER"/{/.}.fasta.gz
