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
