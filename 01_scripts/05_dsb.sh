#!/bin/bash
# Find Double Stranded Breaks (DSBs)

# Parameters
NCPUS=$1

# Global variables
DEDUPLICATED_FOLDER="08_deduplicated"
DSB_FOLDER="09_sites"

# Parallelize on all raw data files
ls -1 -S "$DEDUPLICATED_FOLDER"/*.dedup |
    sort -V |
    parallel -k -j "$NCPUS" ./01_scripts/util/dsb_one_file.py {} "$DSB_FOLDER"/{/.}.dsb
