#!/bin/bash

#SBATCH --job-name=spatial_integration
#SBATCH --output=logs_spatial/glue_%A_%a.out
#SBATCH --error=logs_spatial/glue_%A_%a.err
#SBATCH --array=1-15
#SBATCH --time=1:00:00
#SBATCH --mem=15G
#SBATCH --partition=gpu
#SBATCH --gpus=1


module load conda
conda activate scglue
module load bedtools
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

echo "Running GLUE..."
python scripts/methods/spatial/glue.py $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX

