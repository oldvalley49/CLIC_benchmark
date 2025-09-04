#!/bin/bash

#SBATCH --job-name=r_integration
#SBATCH --output=logs/r_method_%A_%a.out
#SBATCH --error=logs/r_method_%A_%a.err
#SBATCH --array=1-42
#SBATCH --time=5:00:00
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1

module load R

# Read parameters
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/benchmark/parameters/parameters_paired.txt)

# Parse parameters
read DATASET VAR_NUM TOTAL_NUM SUB_NUM COR_METHOD INDEX <<< "$PARAMS"

echo "Running Integration..."
echo "Parameters:"
echo "Dataset: $DATASET"
echo "Variable Features: $VAR_NUM"
echo "Total Number of Features: $TOTAL_NUM"
echo "Subsampling Ratio: $SUB_NUM"
echo "Feature Selection Method: $COR_METHOD"
echo "Replicate Index: $INDEX"

# Define the output directory
OUTDIR="output/paired/rna_atac/${DATASET}/coembed"

# Run Seurat
Rscript --no-save --no-restore scripts/methods/paired/seurat.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX

# Run LIGER
Rscript --no-save --no-restore scripts/methods/paired/liger.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX

# Run BindSC
Rscript --no-save --no-restore scripts/methods/paired/bindsc.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX
rm -r "bindsc-${DATASET}-${VAR_NUM}-${COR_METHOD}-${INDEX}"

echo "All integration methods completed"
