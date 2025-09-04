library(rliger)
library(Seurat)
library(dplyr)
library(ggplot2)
library(class)
library(Matrix)
library(cowplot)
source("scripts/methods/utility.R")
#runLiger_unpaired <- function(rna_batch, atac_batch, rna_dir, atac_dir, var_num, total_num, sub_num, cor_method, index){

args <- commandArgs(trailingOnly = TRUE)
rna_batch <- args[1]
atac_batch <- args[2]
var_num <- as.numeric(args[3])
total_num <- as.numeric(args[4])
cor_method <- args[5]
index <- as.numeric(args[6])
rna_dir <- paste0("data/processed_data/rna_atac/", rna_batch)
atac_dir <- paste0("data/processed_data/rna_atac/", atac_batch)

out_dir <- file.path("output/unpaired/rna_atac", paste0(rna_batch, "+", atac_batch))
plot_dir <- file.path("plots/unpaired/rna_atac", paste0(rna_batch, "+", atac_batch))
dir.create(out_dir, showWarnings=FALSE)
dir.create(plot_dir, showWarnings=FALSE)

# read rna counts matrix as sparse matrix
rna_counts_fp <- file.path(rna_dir, "rna_counts.mtx")
rna_counts <- readMM(rna_counts_fp)
# add cell barcode and gene name to rna_counts
barcodes_rna <- read.csv(file.path(rna_dir, "barcode.csv"))
barcodes_rna <- barcodes_rna$x
genes <- read.csv(file.path(rna_dir, "genes.csv"))
genes <- genes$x
rownames(rna_counts) <- genes
colnames(rna_counts) <- barcodes_rna

# read gene activity score 
activity_counts_fp <- file.path(atac_dir, "activity_counts.mtx")
activity_counts <- readMM(activity_counts_fp)
barcodes_atac <- read.csv(file.path(atac_dir, "barcode.csv"))
barcodes_atac <- barcodes_atac$x
# add cell barcode and gene name to activity counts
activity_genes <- read.csv(file.path(atac_dir, "activity_genes.csv"))
activity_genes <- activity_genes$x
rownames(activity_counts) <- activity_genes
colnames(activity_counts) <- barcodes_atac

# data frame containing correspondance from cell barcode to cell type for RNA
ct_truth_rna <- read.csv(file.path(rna_dir, "ct_annotation.csv"))
ct_truth_rna <- ct_truth_rna$x
barcode_to_annotation_rna <- data.frame(
    barcodes = barcodes_rna,
    ct_truth = ct_truth_rna
)
# data frame containing correspondance from cell barcode to cell type for ATAC
ct_truth_atac <- read.csv(file.path(atac_dir, "ct_annotation.csv"))
ct_truth_atac <- ct_truth_atac$x
barcode_to_annotation_atac <- data.frame(
    barcodes = barcodes_atac,
    ct_truth = ct_truth_atac
)


# initialize seurat object
rna <- CreateSeuratObject(rna_counts, assay="RNA")
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures=var_num)

# filter genes
variable.features <- VariableFeatures(rna)
result <- feature_selection(variable.features, rna_batch, cor_method, var_num, total_num)
use.features <- result$use.features
cor.num <- result$cor.num

# generate run id
run.id <-  paste("liger", rna_batch, atac_batch, cor_method, var_num, cor.num, length(Cells(rna)), ncol(activity_counts), index, sep="_")
run.id

# filter according to the genes

# we won't be able to use genes that don't have fragments in the promoter region and genomic region
use.features = intersect(rownames(activity_counts), use.features)
rna_counts <- rna_counts[use.features, ]
activity_counts <- activity_counts[use.features, ]
# create liger object
liger.obj <- createLiger(list(rna = rna_counts, atac= activity_counts), removeMissing = FALSE)
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
dir.create(file.path(plot_dir, "umap"), showWarnings=FALSE)
ggsave(file.path(plot_dir, "umap",paste0(run.id, ".jpeg")))
# output coembed
coembed <- liger.obj@H.norm

dir.create(file.path(out_dir, "coembed"), showWarnings=FALSE)
write.csv(coembed, file=file.path(out_dir, "coembed", paste0(run.id, ".csv")), row.names = TRUE)

# transfer annotation with KNN classifier(k=30)
rna_cell_num <- ncol(rna)
rna_embed <- coembed[1:rna_cell_num, ]
atac_embed <- coembed[(rna_cell_num+1):nrow(coembed), ]
train_barcodes <- rownames(rna_embed)
train_annotations <- barcode_to_annotation_rna
train_annotations$barcodes <- paste0("rna_", train_annotations$barcodes)
all(train_annotations$barcodes == train_barcodes)
ct.predictions <- knn(train = rna_embed, test = atac_embed, cl = train_annotations$ct_truth, k = 30)
# annotation accuracy table

predictions <- table(barcode_to_annotation_atac$ct_truth, ct.predictions)
predictions <- predictions/rowSums(predictions)  
predictions <- as.data.frame(predictions)

p1 <- ggplot(predictions, aes(Var1, ct.predictions, fill = Freq)) + 
    geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") +
    ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(paste("Variable=", var_num, "Method=", cor_method, "Cell Num=", length(Cells(rna))))

dir.create(file.path(plot_dir, "annotation_accuracy"), showWarnings=FALSE)
dir.create(file.path(out_dir, "ct_annotation"), showWarnings=FALSE)
ggsave(file.path(plot_dir, "annotation_accuracy", paste0(run.id, ".jpeg")))
# save annotation result
cell_annotation_results <- data.frame(barcodes=rownames(atac_embed), predicted.id = ct.predictions, ct_truth=barcode_to_annotation_atac$ct_truth)
write.csv(cell_annotation_results, file.path(out_dir,"ct_annotation", paste0(run.id, ".csv")))

   