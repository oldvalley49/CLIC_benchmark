#!/bin/bash

#SBATCH --job-name=r_integration
#SBATCH --output=logs/paired/paired_cpu_%A_%a.out
#SBATCH --error=logs/paired/paired_cpu_%A_%a.err
#SBATCH --array=1-672
#SBATCH --time=4:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1

module load R

# Read parameters
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/benchmark/parameters/parameters_paired.txt)

# Parse parameters
read DATASET VAR_NUM TOTAL_NUM SUB_NUM COR_METHOD INDEX ACTIVITY_MODEL <<< "$PARAMS"

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
Rscript --no-save --no-restore scripts/methods/paired/seurat.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX $ACTIVITY_MODEL

# Run LIGER
Rscript --no-save --no-restore scripts/methods/paired/liger.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX $ACTIVITY_MODEL

# Run BindSC
Rscript --no-save --no-restore scripts/methods/paired/bindsc.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX $ACTIVITY_MODEL
rm -r "bindsc-${DATASET}-${VAR_NUM}-${COR_METHOD}-${ACTIVITY_MODEL}-${INDEX}"

echo "All integration methods completed"
