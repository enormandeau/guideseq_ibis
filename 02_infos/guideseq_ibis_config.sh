#!/bin/bash

# Modify the following parameter values according to your experiment
# Do not modify the parameter names or remove parameters
# Do not add spaces around the equal (=) sign
# It is a good idea to try to run Barque with different parameters 

# Global parameters
NCPUS=40                    # Number of CPUs to use. A lot of the steps are parallelized (int, 1+)

MIN_LENGTH=100              # Length of fragments used for DBS (int, 50+, 100 works well)
CROP_LENGTH=200             # Cut reads to this length after filtering (int, => MIN_LENGTH)
GENOME_FILE="genome.fasta"  # Name of genome file in 03_genome folder

# Double stranded breaks (DSB) parameters
MIN_COVERAGE=20             # Minimum total coverage to detect a target or off-target
                            #   Recommendation is around 10 * (depth of coverage / 100,000)
                            #   However, this will vary as a function of the poportion of target and
                            #   off-target hits compared to background noise, as well as desired sensitivity.
WINDOW_SIZE=10              # Window size on each side of highest point of coverage to count alignments
POSITION_ERROR=5            # Accept DSBs with forward and reverse targets at slightly different positions
BIN_SIZE=10000              # Divide chromosomes in bins of that size to locate genes. Keep at 1000 or 10000
