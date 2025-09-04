#!/bin/bash

#SBATCH --job-name=benchmark
#SBATCH --output=logs/benchmark_%A_%a.out
#SBATCH --error=logs/benchmark_%A_%a.err
#SBATCH --array=1-7
#SBATCH --time=5:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

module load conda
conda activate scglue

PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/benchmark/paired_tissues.txt)

read TISSUE <<< "$PARAMS"

echo "Running Benchmark..."
echo "Tissue: $TISSUE"

python scripts/benchmark/run_benchmark.py $TISSUE