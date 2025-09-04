library(Seurat)
library(Matrix)
library(glue)
source("scripts/preprocessing/utility.R")

tissues <- c("mBrain2")


for (tissue in tissues){
    data_dir <- paste0("data/processed_data/spatial/", tissue)
    rna_counts <- load_rna(data_dir)
    rna <- CreateSeuratObject(rna_counts, assay="RNA")
    rna$orig.ident <- "RNA"
    rna <- NormalizeData(rna)
    refdata <- GetAssayData(rna, assay = "RNA", layer = "data")
    fp <- glue("data/processed_data/spatial/{tissue}/norm_rna_counts.RDS")
    saveRDS(refdata, file=fp)
}
