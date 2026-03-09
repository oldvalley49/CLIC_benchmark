library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(glue)

output_dir <- "output/spatial"
tissues <- list.dirs(output_dir, full.names = FALSE, recursive = FALSE)

for (tissue in tissues) {
    if (tissue == "mBrain2") {
        marker_genes <- c("Cd63")
    } else {
        next
    }
    
    gex_dir <- glue("{output_dir}/{tissue}/geximpute")
    methods <- list.files(gex_dir)
    
    for (method in methods) {
        print(method)
        method_name <- strsplit(method, "\\.")[[1]][1]
        
        gex <- readRDS(glue("{gex_dir}/{method}"))
        colnames(gex) <- substr(colnames(gex), 1, 16)
        
        seurat.obj <- CreateSeuratObject(gex, assay = "SPATIAL")
        seurat.obj <- NormalizeData(seurat.obj)
        seurat.obj[["SPATIAL"]]$data <- gex
        
        image <- Read10X_Image(image.dir = glue("data/processed_data/spatial/{tissue}/spaceranger"), filter.matrix = TRUE)
        DefaultAssay(image) <- "SPATIAL"
        seurat.obj[["slice1_"]] <- image
        
        for (gene in marker_genes) {
            p <- SpatialFeaturePlot(seurat.obj, features = gene, pt.size.factor = 4)
            
            ggsave(filename = glue("plots/spatial/{tissue}/spatial/{method_name}_{gene}.svg"), plot = p, create.dir=TRUE)
        }
    }
}
