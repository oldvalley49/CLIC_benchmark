#!/bin/bash

#SBATCH --job-name=spatial_corr
#SBATCH --output=logs/spatial_corr%A_%a.out
#SBATCH --error=logs/spatial_corr_%A_%a.err
#SBATCH --array=1-15
#SBATCH --time=2:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1

module load R
TISSUE="mBrain2"
# Read parameters
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" output/spatial/"${TISSUE}"/gex_impute_methods.txt)

# Parse parameters
read METHOD <<< "$PARAMS"

echo "Computing Correlation"
echo "Method: $METHOD"

Rscript --no-save --no-restore scripts/benchmark/imputation_correlation_spatial.R $TISSUE $METHOD

echo "All integration methods completed"
