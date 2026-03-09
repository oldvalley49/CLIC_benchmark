library(rliger)
library(Seurat)
library(dplyr)
library(ggplot2)
library(class)
library(Matrix)
library(cowplot)
library(FNN)
library(devtools)
devtools::load_all('/dcs07/hongkai/data/tomo/CLIC')
source("scripts/methods/utility.R")

# parse input argument
args <- commandArgs(trailingOnly = TRUE)
tissue = args[1]
var_num <- as.numeric(args[2])
total_num <- as.numeric(args[3])
sub_num <- as.numeric(args[4])
cor_method <- args[5]
index <- as.numeric(args[6])
activity_model <- args[7]
data_dir <- paste0("data/processed_data/rna_atac/", tissue)
dir.create(file.path("output/paired/rna_atac", tissue), showWarnings=FALSE)
dir.create(file.path("plots/paired/rna_atac", tissue), showWarnings=FALSE)
# load data
rna_counts <- load_rna(data_dir)
atac_counts <- load_atac(data_dir)
if (activity_model == 'signac') {
    activity_counts = load_activity(data_dir)
} else {
    activity_counts = load_maestro(data_dir)
}
barcode_to_annotation <- load_annotations(data_dir)

activity_counts <- activity_counts[intersect(rownames(activity_counts), 
                                             rownames(rna_counts)), ]
# subsample
subsampled <- subsample(rna_counts, atac_counts, activity_counts, barcode_to_annotation, sub_num)
rna_counts <- subsampled$rna_counts
atac_counts <- subsampled$atac_counts
activity_counts <- subsampled$activity_counts
barcode_to_annotation <- subsampled$barcode_to_annotation

rm(atac_counts)
# initialize seurat object
rna <- CreateSeuratObject(rna_counts, assay="RNA")
rna <- NormalizeData(rna)

if (startsWith(tissue, "m")) {
    out <- FindCLICFeatures(rna, score_name = paste0('mouse-', cor_method), initial_variable_features_num = var_num, nfeatures = total_num)
} else {
    out <- FindCLICFeatures(rna, score_name = paste0('human-', cor_method), initial_variable_features_num = var_num, nfeatures = total_num)
} 

use.features <- out$use_features
cor.num <- out$clic_num

# generate run id
run.id <-  paste("liger", cor_method, var_num, cor.num, length(Cells(rna)), activity_model, index, sep="_")
print(run.id)

# filter according to the genes

# we won't be able to use genes that don't have fragments in the promoter region and genomic region

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
p1 <- plotDatasetDimRed(liger.obj, title=paste("Variable=", var_num, "Method=", cor_method, "Cell Num=", length(Cells(rna))))
dir.create(file.path("plots/paired/rna_atac", tissue, "umap"), showWarnings=FALSE)
ggsave(file.path("plots/paired/rna_atac", tissue, "umap",paste0(run.id, ".jpeg")))
# output coembed
coembed <- liger.obj@H.norm
dir.create(file.path("output/paired/rna_atac", tissue, "coembed"))
write.csv(coembed, file=file.path("output/paired/rna_atac", tissue, "coembed", paste0(run.id, ".csv")), row.names = TRUE)

# transfer annotation with KNN classifier(k=30)
cell_num <- ncol(rna)
rna_embed <- coembed[1:cell_num, ]
atac_embed <- coembed[(cell_num+1):nrow(coembed), ]
train_barcodes <- rownames(rna_embed)
train_annotations <- barcode_to_annotation
train_annotations$barcodes <- paste0("rna_", train_annotations$barcodes)
all(train_annotations$barcodes == train_barcodes)
ct.predictions <- knn(train = rna_embed, test = atac_embed, cl = train_annotations$ct_truth, k = 30)

# annotation accuracy table
predictions <- table(train_annotations$ct_truth, ct.predictions)
predictions <- predictions/rowSums(predictions)  
predictions <- as.data.frame(predictions)

p1 <- ggplot(predictions, aes(Var1, ct.predictions, fill = Freq)) + 
    geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") +
    ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(paste("Variable=", var_num, "Method=", cor_method, "Cell Num=", length(Cells(rna))))
dir.create(file.path("plots/paired/rna_atac", tissue, "annotation_accuracy"))
ggsave(file.path("plots/paired/rna_atac", tissue, "annotation_accuracy", paste0(run.id, ".jpeg")))
# save annotation result
cell_annotation_results <- data.frame(barcodes=rownames(atac_embed), predicted.id = ct.predictions, ct_truth=train_annotations$ct_truth)
dir.create(file.path("output/paired/rna_atac", tissue, "ct_annotation"), showWarnings=FALSE)
write.csv(cell_annotation_results, file.path("output/paired/rna_atac", tissue, "ct_annotation", paste0(run.id, ".csv")))
