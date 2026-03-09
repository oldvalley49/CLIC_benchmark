library(rliger)
library(Seurat)
library(dplyr)
library(ggplot2)
library(class)
library(Matrix)
library(cowplot)
library(FNN)
source("scripts/methods/utility.R")

# parse input argument
args <- commandArgs(trailingOnly = TRUE)
tissue = args[1]
var_num <- as.numeric(args[2])
total_num <- as.numeric(args[3])
sub_num <- as.numeric(args[4])
cor_method <- args[5]
index <- as.numeric(args[6])
data_dir <- paste0("data/processed_data/spatial/", tissue)
dir.create(file.path("output/spatial", tissue), showWarnings=FALSE)
dir.create(file.path("plots/spatial", tissue), showWarnings=FALSE)
# load data
rna_counts <- load_rna(data_dir)
atac_counts <- load_atac(data_dir)
activity_counts <- load_activity(data_dir)

# subsample
subsampled <- subsample_no_annotation(rna_counts, atac_counts, activity_counts, sub_num)
rna_counts <- subsampled$rna_counts
atac_counts <- subsampled$atac_counts
activity_counts <- subsampled$activity_counts

rm(atac_counts)
# initialize seurat object
rna <- CreateSeuratObject(rna_counts, assay="RNA")
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures=var_num)

# filter genes
variable.features <- VariableFeatures(rna)
result <- feature_selection(variable.features, tissue, cor_method, var_num, total_num)
use.features <- result$use.features
cor.num <- result$cor.num
# generate run id
run.id <-  paste("liger", cor_method, var_num, cor.num, length(Cells(rna)), index, sep="_")
run.id

use.features = intersect(rownames(activity_counts), use.features)


# convert to dgC matrix
rna_counts <- as(rna_counts, "CsparseMatrix")
activity_counts <- as(activity_counts, "CsparseMatrix")
# create liger object
liger.obj <- createLiger(list(rna = rna_counts, atac = activity_counts), removeMissing=TRUE)
liger.obj <- normalize(liger.obj)
#liger.obj <- selectGenes(liger.obj, useDatasets = "rna")
liger.obj@varFeatures <- use.features
liger.obj <- scaleNotCenter(liger.obj)
liger.obj <- runIntegration(liger.obj, k = 20)
liger.obj <- quantileNorm(liger.obj)
liger.obj <- runCluster(liger.obj, nNeighbors = 30, resolution = 0.2)
liger.obj <- runUMAP(liger.obj, nNeighbors = 30, minDist = 0.3)

# plot UMAP
p1 <- plotDatasetDimRed(liger.obj, title=paste("Variable=", var_num, "Method=", cor_method, "Cell Num=", length(colnames(rna_counts))))
dir.create(file.path("plots/spatial", tissue, "umap"), showWarnings=FALSE)
ggsave(file.path("plots/spatial", tissue, "umap",paste0(run.id, ".jpeg")))
# output coembed
coembed <- liger.obj@H.norm
dir.create(file.path("output/spatial", tissue, "coembed"))
write.csv(coembed, file=file.path("output/spatial", tissue, "coembed", paste0(run.id, ".csv")), row.names = TRUE)
