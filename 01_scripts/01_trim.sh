#!/bin/bash
# Trim reads and remove
# - Reads with bad quality
# - Short reads

# Parameters
NCPUS=$1
MIN_LENGTH=$2
CROP_LENGTH=$3

# Global variables
DATA_FOLDER="04_data"

# Parallelize on all raw data files
ls -1 -S "$DATA_FOLDER"/*_R1_001.fastq.gz |
    perl -pe 's/_R[12]_001\.fastq\.gz/_/' |
    parallel -k -j "$NCPUS" ./01_scripts/util/trimmomatic.sh {} "$MIN_LENGTH" "$CROP_LENGTH"
