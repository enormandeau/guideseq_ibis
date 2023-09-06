#!/bin/bash

# Modify the following parameter values according to your experiment
# Do not modify the parameter names or remove parameters
# Do not add spaces around the equal (=) sign
# It is a good idea to try to run Barque with different parameters 

# Global parameters
NCPUS=40                    # Number of CPUs to use. A lot of the steps are parallelized (int, 1+)

MIN_LENGTH=100              # Length to trim sequences (int, 50+)
CROP_LENGTH=200             # Cut reads to this length after filtering (int, 100+)

GENOME_FILE="genome.fasta"  # Name of genome file in 03_genome folder

# Double stranded breaks (DSB) parameters
MIN_COVERAGE=20             # Minimum total coverage to detect a target or off-target
WINDOW_SIZE=10              # Window size on each side of highest point of coverage to count alignments
POSITION_ERROR=5            # Accept DSBs with forward and reverse targets at slightly different positions
BIN_SIZE=10000              # Divide chromosomes in bins of that size to locate genes. Keep at 1000 or 10000
