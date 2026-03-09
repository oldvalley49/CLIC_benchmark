#!/bin/bash

#SBATCH --job-name=r_unpaired_integration
#SBATCH --output=logs_unpaired/r_method_%A_%a.out
#SBATCH --error=logs_unpaired/r_method_%A_%a.err
#SBATCH --array=11-45
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

module load R

# Read parameters
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/benchmark/parameters/parameters_unpaired.txt)

# Parse parameters
read RNA_BATCH ATAC_BATCH VAR_NUM TOTAL_NUM COR_METHOD INDEX <<< "$PARAMS"

echo "Running Integration..."
echo "Parameters:"
echo "GEX Batch: $RNA_BATCH"
echo "ATAC Batch: $ATAC_BATCH"
echo "Number of Variable Features: $VAR_NUM"
echo "Feature Selection Method: $COR_METHOD"
echo "Replicate Index: $INDEX"

#Run Seurat
echo "Running Seurat..."
Rscript --no-save --no-restore scripts/methods/unpaired/seurat.R $RNA_BATCH $ATAC_BATCH $VAR_NUM $TOTAL_NUM $COR_METHOD $INDEX

# #Run LIGER
# echo "Running LIGER..."
# Rscript --no-save --no-restore scripts/methods/unpaired/liger.R $RNA_BATCH $ATAC_BATCH $VAR_NUM $TOTAL_NUM $COR_METHOD $INDEX

# #Run bindSC
# echo "Running bindSC..."
# Rscript --no-save --no-restore scripts/methods/unpaired/bindsc.R $RNA_BATCH $ATAC_BATCH $VAR_NUM $TOTAL_NUM $COR_METHOD $INDEX
# rm -r "bindsc-${RNA_BATCH}-${ATAC_BATCH}-${VAR_NUM}-${COR_METHOD}-${INDEX}"

echo "All integration methods completed"
