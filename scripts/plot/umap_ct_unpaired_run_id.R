library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(cowplot)
library(glue)
library(stringr)
source("scripts/utility/infer_barcode.R")
### INFO
rna_tissue <- "BMMC-s1d2"
atac_tissue <- "BMMC-s2d4"
run.id <- "liger_BMMC-s1d2_BMMC-s2d4_corr_5000_6533_6740_6111_1.csv"

### SCRIPT
combined <- glue("{rna_tissue}+{atac_tissue}")
dir.create(glue("plots/unpaired/rna_atac/{combined}/umap_ct"), showWarnings=FALSE)
tested_combination <- strsplit(run.id, "\\.")[[1]][1]

# load coembed
coembed <- read.csv(glue("output/unpaired/rna_atac/{combined}/coembed/{run.id}"), row.names = 1, header=TRUE)
coembed <- as.matrix(coembed)

idx <- sub(".*?(\\d+)$", "\\1", colnames(coembed))
colnames(coembed) <- paste0("PC_", idx)


# load cell-type labels

ct_rna <- read.csv(glue("data/processed_data/rna_atac/{rna_tissue}/ct_annotation.csv"), row.names=1)
ct_atac <- read.csv(glue('data/processed_data/rna_atac/{atac_tissue}/ct_annotation.csv'), row.names=1)

ct_df <- rbind(ct_rna, ct_atac)
colnames(ct_df) <- "celltype"
rownames(coembed) <- sapply(rownames(coembed), infer_barcode)


# DUMMY count matrix
dummy_counts <- Matrix::Matrix(data = 0, nrow = 2, ncol = nrow(coembed), sparse = TRUE)
rownames(dummy_counts) <- c("FAKE1", "FAKE2")
colnames(dummy_counts) <- rownames(coembed)
seurat_obj <- CreateSeuratObject(counts = dummy_counts, meta.data = ct_df)
pca_reduction <- CreateDimReducObject(embeddings = coembed, assay = "RNA", key="PC_")
seurat_obj[["pca"]] <- pca_reduction
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", seed.use = 42,dims = 1:ncol(coembed))
p <- DimPlot(seurat_obj, reduction = "umap", group.by = "celltype", shuffle=TRUE) + coord_fixed()
ggsave(
    filename = glue("plots/unpaired/rna_atac/{combined}/umap_ct/{run.id}.svg"),
    plot = p,
    width = 10,
    height = 10
    )

