#!/bin/bash
#
#SBATCH --mem=30G
#SBATCH --time=5:00:00

module load R
module load samtools
module load htslib

echo "`date`: Running R Job"
R CMD BATCH scripts/preprocessing/mBrain_preprocess.R

echo "`date`: R job complete"
exit 0