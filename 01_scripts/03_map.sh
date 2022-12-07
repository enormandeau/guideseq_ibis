#!/bin/bash
# Map fasta reads to genome

# Parameters
NCPUS=$1
GENOME_FILE="$2"

# Global variables
GENOME_FOLDER="03_genome"
EXTRACTED_FOLDER="06_extracted"

# Parallelize on all data files
ls -1 -S "$EXTRACTED_FOLDER"/*.fasta.gz |
    perl -pe 's/\.fasta\.gz//' |
    parallel -k -j 1 ./01_scripts/util/map_one_file.sh {} "$NCPUS" "$GENOME_FILE"
