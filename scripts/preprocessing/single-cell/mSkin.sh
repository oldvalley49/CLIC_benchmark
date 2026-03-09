#!/bin/bash
#
#SBATCH --output=logs/preprocessing/mSkin_%A_%a.out
#SBATCH --error=logs/preprocessing/mSkin_%A_%a.err
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

module load R
module load samtools
module load htslib

Rscript --no-save --no-restore scripts/preprocessing/single-cell/mSkin.R
