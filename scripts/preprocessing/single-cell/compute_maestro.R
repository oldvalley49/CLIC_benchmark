library(Seurat)
library(Matrix)
library(glue)
source("scripts/preprocessing/utility.R")

human_tissues <- c("BMMC-s1d2", "BMMC-s4d8", "PBMC", "TDBM")
mouse_tissues <- c("mSkin") # "mBrain", "mRetina", 
 
for (tissue in human_tissues){
    data_dir <- paste0("data/processed_data/rna_atac/", tissue)
    compute_MAESTRO(data_dir, 'GRCh38')
}

for (tissue in mouse_tissues){
    data_dir <- paste0("data/processed_data/rna_atac/", tissue)
    compute_MAESTRO(data_dir, 'GRCm38')
}
