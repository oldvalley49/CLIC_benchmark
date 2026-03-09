#!/bin/bash

#SBATCH --job-name=glue_integration
#SBATCH --output=logs/glue_%A_%a.out
#SBATCH --error=logs/glue_%A_%a.err
#SBATCH --array=10-336
#SBATCH --time=4:00:00
#SBATCH --mem=70G
#SBATCH --partition=gpu
#SBATCH --gpus=1

module load conda
module load bedtools
conda activate scglue


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

OUTDIR="output/paired/rna_atac/${DATASET}/coembed"

# Run GLUE only if output doesn't already exist
echo "Running GLUE..."
python scripts/methods/paired/glue.py $DATASET $VAR_NUM $TOTAL_NUM $SUB_NUM $COR_METHOD $INDEX $ACTIVITY_MODEL

