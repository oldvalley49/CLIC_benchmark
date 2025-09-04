# Benchmark Pipeline for CLIC scores

## Overview
The CLIC package provides the tool for feature selection for integrating unpaired scRNA-seq and scATAC-seq data in a way that leverages multiome data from ENCODE. 

CLIC package can be found [here](https://github.com/oldvalley49/CLIC)

This repository contains all of the scripts used to evaluate CLIC against baseline variable feature selection when used within the workflow of SOTA methods. 

For more details about the algorithm as well as the benchmarking process, please refer to our preprint: PREPRINT LINK PENDING

## Description

`scripts` directory contains all relevant scripts used to produce benchmark results. 
`scripts/preprocessing`: preprocessing scripts used to extract data from raw data files downloaded from original publications. 
`scripts/methods`: scripts for running each integration algorithm (Seurat, GLUE, Liger, bindSC) with CLIC/baseline feature selection
`scripts/benchmark`: SLURM job submission scripts for running benchmark experiments on HPC. In particular:
`scripts/benchmark/submit_{paired, unpaired, spatial}.sh` and `scripts/benchmark/submit_gpu_{paired, unpaired, spatial}` contain submission scripts for running integration methods. 
`scripts/benchmark/parameters` contain code used to generate parameter pairings used to assess robustness to parameter settings
`scripts/run_benchmark_{'', 'unpaired', 'spatial'}` contain scripts used to compute key metrics used to assess integration quality. 
