#!/bin/bash
#
#SBATCH --output=logs/preprocessing/mBrain_%A_%a.out
#SBATCH --error=logs/preprocessing/mBrain_%A_%a.err
#SBATCH --mem=30G
#SBATCH --time=5:00:00

module load R
module load samtools
module load htslib

Rscript --no-save --no-restore scripts/preprocessing/single-cell/mBrain.R
