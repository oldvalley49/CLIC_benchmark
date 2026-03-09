#!/bin/bash

#SBATCH --job-name=glue_unpaired_integration
#SBATCH --output=logs_unpaired/glue_%A_%a.out
#SBATCH --error=logs_unpaired/glue_%A_%a.err
#SBATCH --array=1-15
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --partition=gpu
#SBATCH --gpus=1

module load conda
module load bedtools
conda activate scglue

# Read parameters
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/benchmark/parameters/parameters_unpaired.txt)

# Parse parameters
read RNA_BATCH ATAC_BATCH VAR_NUM TOTAL_NUM COR_METHOD INDEX <<< "$PARAMS"

echo "Running Integration..."
echo "Parameters:"
echo "GEX Batch: $RNA_BATCH"
echo "ATAC Batch: $ATAC_BATCH"
echo "Number of Variable Features: $VAR_NUM"
echo "Total Number of Features: $TOTAL_NUM"
echo "Feature Selection Method: $COR_METHOD"
echo "Replicate Index: $INDEX"


echo "Running GLUE..."
python scripts/methods/unpaired/glue.py $RNA_BATCH $ATAC_BATCH $VAR_NUM $TOTAL_NUM $COR_METHOD $INDEX