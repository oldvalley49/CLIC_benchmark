#!/bin/bash
#
#SBATCH --mem=30G
#SBATCH --time=4:00:00

module load conda
conda activate scglue

python scripts/benchmark/run_benchmark_spatial.py