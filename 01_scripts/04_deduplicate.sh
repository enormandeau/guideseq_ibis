#!/bin/bash
# Remove PCR duplicates

# Parameters
NCPUS=$1

# Global variables
ALIGNED_FOLDER="07_aligned" 
DEDUPLICATED_FOLDER="08_deduplicated"

# Parallelize on all raw data files
ls -1 -S "$ALIGNED_FOLDER"/*.sam |
    parallel -k -j "$NCPUS" ./01_scripts/util/deduplicate_one_file.py {} "$DEDUPLICATED_FOLDER"/{/.}.dedup
