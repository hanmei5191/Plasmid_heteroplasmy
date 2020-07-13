#!/bin/bash

#SBATCH --job-name=Plasmid_chr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hzm5191@psu.edu
set -x

Rscript optimized.R
