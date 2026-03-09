#!/bin/bash

#SBATCH --job-name=paired_benchmark
#SBATCH --output=logs/benchmark/paired_%A_%a.out
#SBATCH --error=logs/benchmark/paired_%A_%a.err
#SBATCH --array=7
#SBATCH --time=2-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1

module load conda
conda activate scglue

PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" scripts/benchmark/compute_metrics/paired_tissues.txt)

read TISSUE <<< "$PARAMS"

echo "Running Benchmark..."
echo "Tissue: $TISSUE"

python scripts/benchmark/compute_metrics/run_benchmark.py $TISSUE