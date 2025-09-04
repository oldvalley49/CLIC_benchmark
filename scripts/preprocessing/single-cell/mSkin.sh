#!/bin/bash
#
#SBATCH --mem=100G
#SBATCH --time=2-00:00:00

module load R
module load samtools
module load htslib

echo "`date`: Running R Job"
R CMD BATCH scripts/preprocessing/mSkin_preprocess.R

echo "`date`: R job complete"
exit 0