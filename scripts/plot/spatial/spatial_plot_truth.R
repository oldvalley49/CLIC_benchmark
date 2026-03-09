library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(glue)

data_dir <- "data/processed_data/spatial"
tissues <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

for (tissue in tissues) {
    if (tissue == "mBrain2") {
        marker_genes <- c("Bcl11b")
    } else {
        next
    }
    tissue_dir <- glue("{data_dir}/{tissue}")
    gex <- readRDS(glue("{tissue_dir}/norm_rna_counts.RDS"))
    colnames(gex) <- substr(colnames(gex), 1, 16)
    
    seurat.obj <- CreateSeuratObject(gex, assay = "SPATIAL")
    seurat.obj <- NormalizeData(seurat.obj)
    seurat.obj[["SPATIAL"]]$data <- gex
    
    image <- Read10X_Image(image.dir = glue("{tissue_dir}/spaceranger"), filter.matrix = TRUE)
    DefaultAssay(image) <- "SPATIAL"
    seurat.obj[["slice1_"]] <- image
    
    for (gene in marker_genes) {
        p <- SpatialFeaturePlot(seurat.obj, features = gene, pt.size.factor = 4)
        ggsave(filename = glue("{tissue_dir}/spatial_plot_{gene}.svg"), plot = p)
    }
}
