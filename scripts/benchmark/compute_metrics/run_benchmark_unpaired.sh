#!/bin/bash
#
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=logs/bench_UNP_%A_%a.out

module load conda
conda activate scglue

python scripts/benchmark/run_benchmark_unpaired.py