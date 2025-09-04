library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(cowplot)
library(glue)
library(stringr)

### INFO
tissue <- "BMMC-s4d8"
run.id <- "liger_corr_5000_7599_9876_1.csv"

### SCRIPT
dir.create(glue("plots/paired/rna_atac/{tissue}/umap_ct"), showWarnings=FALSE)
tested_combination <- strsplit(run.id, "\\.")[[1]][1]

# load coembed
coembed <- read.csv(glue("output/paired/rna_atac/{tissue}/coembed/{run.id}"), header=TRUE)
coembed[,1] <- make.unique(coembed[,1])
rownames(coembed) <- coembed[,1]
coembed <- coembed[, -1]
coembed <- as.matrix(coembed)

idx <- sub(".*?(\\d+)$", "\\1", colnames(coembed))
colnames(coembed) <- paste0("PC_", idx)

cell_num <- strsplit(run.id, "_")[[1]][5]
# load cell-type labels
ct <- read.csv(glue("output/paired/rna_atac/{tissue}/ct_annotation/{run.id}"), row.names=1)
ct <- c(ct$ct_truth, ct$ct_truth)
# fix duplicate row names
names(ct) <- rownames(coembed)

ct_df <- data.frame(celltype = ct)
rownames(ct_df) <- names(ct)
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
    filename = glue("plots/paired/rna_atac/{tissue}/umap_ct/{tested_combination}.svg"),
    plot = p,
    device = 'svg',
    )

