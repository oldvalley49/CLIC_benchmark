#!/bin/bash

#SBATCH --job-name=spatial_integration
#SBATCH --output=logs_spatial/r_method_%A_%a.out
#SBATCH --error=logs_spatial/r_method_%A_%a.err
#SBATCH --array=1-15
#SBATCH --time=2:00:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1

module load R

# Read parameters
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/benchmark/parameters/parameters_paired_spatial.txt)

# Parse parameters
read DATASET VAR_NUM TOTAL_NUM SUB_NUM COR_METHOD INDEX <<< "$PARAMS"

echo "Running Spatial Integration..."
echo "Parameters:"
echo "Dataset: $DATASET"
echo "Variable Features: $VAR_NUM"
echo "Total Number of Features: $TOTAL_NUM"
echo "Subsampling Ratio: $SUB_NUM"
echo "Correlation Method: $COR_METHOD"
echo "Replicate Index: $INDEX"

# Run Seurat
echo "Running Seurat..."
Rscript --no-save --no-restore scripts/methods/spatial/seurat.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX

# #Run LIGER
echo "Running LIGER..."
Rscript --no-save --no-restore scripts/methods/spatial/liger.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX

# Run bindSC
echo "Running bindSC..."
Rscript --no-save --no-restore scripts/methods/spatial/bindsc.R $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX
rm -r "bindsc-${DATASET}-${VAR_NUM}-${COR_METHOD}-${INDEX}"

echo "All integration methods completed"
