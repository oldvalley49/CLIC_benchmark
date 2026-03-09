# Benchmark Pipeline for CLIC scores

## Overview
The CLIC package provides the tool for feature selection for integrating unpaired scRNA-seq and scATAC-seq data in a way that leverages multiome data from ENCODE. 

CLIC package can be found [here](https://github.com/oldvalley49/CLIC)

This repository contains all of the scripts used to evaluate CLIC against baseline variable feature selection when used within the workflow of SOTA methods. 

For more details about the algorithm as well as the benchmarking process, please refer to our preprint: 

PREPRINT LINK PENDING

## Description

The `scripts` directory contains all code used to generate benchmark results:  

- **`scripts/preprocessing/`**  
  Scripts for preprocessing raw data downloaded from original publications.  

- **`scripts/methods/`**  
  Scripts for running each integration algorithm (Seurat, GLUE, Liger, bindSC) with either CLIC-based or baseline variable feature selection.  
  - `paired` includes method scripts for paired scenarios
  - `unpaired` include method scripts for unpaired scenarios
  - `spatial` include method scripts for spatial scenarios
- **`scripts/benchmark/`**  
  Contains SLURM job submission scripts and parameter configurations for running benchmark experiments on HPC:  
  - `submit_jobs/submit_{paired, unpaired, spatial}.sh` and `submit_jobs/submit_gpu_{paired, unpaired, spatial}.sh`: submission scripts for integration methods. Paired, unpaired, and spatial each represent different benchmarking scenarios in the manuscript. 
  - `parameters/`: code for generating parameter combinations used to assess robustness across settings.  Output textfiles are used to submit array jobs with different parameter settings
  - `compute_metrics/run_benchmark_{'', 'unpaired', 'spatial'}`: scripts for computing key metrics to evaluate integration quality.  

- **`scripts/plot/`**  
  Scripts for plotting figures used in publication
