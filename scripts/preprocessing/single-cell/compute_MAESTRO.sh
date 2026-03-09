#!/bin/bash
#
#SBATCH --job-name=computeMAESTRO
#SBATCH --mem=50G
#SBATCH --time=6:00:00
#SBATCH --output=logs/compute_MAESTRO/r_method_%A_%a.out
#SBATCH --error=logs/compute_MAESTRO/r_method_%A_%a.err

module load R

Rscript --no-save --no-restore scripts/preprocessing/single-cell/compute_maestro.R