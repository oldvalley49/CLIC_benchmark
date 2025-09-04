#runbindSC_unpaired <- function(rna_batch, atac_batch, rna_dir, atac_dir, var_num, total_num, sub_num, cor_method, index){
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(bindSC)
library(Matrix)
library(irlba)
library(umap)
library(Signac)
library(class)
library(cowplot)
library(glue)
source("scripts/methods/utility.R")

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

# read atac counts as sparse matrix
atac_counts_fp <- file.path(atac_dir, "atac_counts.mtx")
atac_counts <- readMM(atac_counts_fp)
# add cell barcode and gene name to rna_counts
barcodes_atac <- read.csv(file.path(atac_dir, "barcode.csv"))
barcodes_atac <- barcodes_atac$x
peaks <- read.csv(file.path(atac_dir, "peaks.csv"))
peaks <- peaks$x
rownames(atac_counts) <- peaks
colnames(atac_counts) <- barcodes_atac


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

# perform subsampling
# perform subsampling
# num_total <- length(barcodes)
# num_subsample <- sub_num * num_total%/%4
# barcode_subsampled <- sample(barcodes, num_subsample, replace=FALSE)
# rna_counts <- rna_counts[, barcode_subsampled]
# atac_counts <- atac_counts[, barcode_subsampled]
# activity_counts <- activity_counts[, barcode_subsampled]
# barcode_to_annotation <- barcode_to_annotation %>%
#    dplyr::filter(barcodes %in% barcode_subsampled) %>%
#    mutate(order = match(barcodes, barcode_subsampled)) %>%
#    arrange(order) %>%
#    dplyr::select(-order)

# initialize seurat object RNA data and perform preprocessing including clustering step required for bindsc
rna <- CreateSeuratObject(rna_counts, assay="RNA")
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures=var_num)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20, reduction="pca")
rna <- FindClusters(rna, resolution = 0.5)
colnames(rna) <- paste0("RNA_", colnames(rna))
# initialize seurat object for ATAC data and perform preprocessing including clustering step required for bindsc
colnames(atac_counts) <- paste0("ATAC_", colnames(atac_counts))
colnames(activity_counts) <- paste0("ATAC_", colnames(activity_counts))
# granges.counts <- StringToGRanges(rownames(atac_counts))
# granges.use <- seqnames(granges.counts) %in% standardChromosomes(granges.counts)
# atac_counts <- atac_counts[as.vector(granges.use), ]
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"
chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'hg38',
    # min.cells = 10,
    # annotation = annotations
)
atac <- CreateSeuratObject(counts=chrom_assay, assay="ATAC")
atac$orig.ident <- "ATAC"
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- FindNeighbors(atac, dims = 1:20, reduction = "lsi")
atac <- FindClusters(atac, resolution = 0.5)
# filter genes
variable.features <- VariableFeatures(rna)
result <- feature_selection(variable.features, rna_batch, cor_method, var_num, total_num)
use.features <- result$use.features
cor.num <- result$cor.num
# generate run id
run.id <-  paste("bindsc", rna_batch, atac_batch, cor_method, var_num, cor.num, length(Cells(rna)), length(Cells(atac)), index, sep="_")
run.id

# we won't be able to use genes that don't have fragments in the promoter region and genomic region
atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity_counts)
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)

use.features = intersect(rownames(activity_counts), use.features)
X <- rna[["RNA"]]$data[use.features, ]
Z0 <- atac[["ACTIVITY"]]$data[use.features,]
# Y <- atac_counts
# dimensional reduction
out <- dimReduce( dt1 =  X, dt2 = Z0,  K = 30)
x <- out$dt1
z0 <- out$dt2
y  <- atac@reductions$lsi@cell.embeddings

# run integration
res <- BiCCA( X = t(x) ,
                Y = t(y), 
                Z0 =t(z0), 
                X.clst = rna$seurat_clusters,
                Y.clst = atac$seurat_clusters,
                alpha = 0.5, 
                lambda = 0.5,
                K = 15,
                temp.path  = glue("bindsc-{rna_batch}-{atac_batch}-{var_num}-{cor_method}-{index}"),
                num.iteration = 50,
                tolerance = 0.01,
                save = TRUE,
                parameter.optimize = FALSE,
                block.size = 0)

# plot coembedded space
dim(res$u)
dim(res$r)
umap_plt <- umap(rbind(res$u, res$r))
umap_plt  <- data.frame("UMAP1"=umap_plt$layout[,1],
                        "UMAP2"=umap_plt$layout[,2],
                        "data" = c(rep("scRNA-seq",ncol(X)),
                                    rep("scATAC-seq",ncol(atac_counts))))

p <- ggplot(umap_plt, aes(x = UMAP1, y = UMAP2, color = data)) + 
    geom_point(size = 0.1) +
    theme_minimal() +
    ggtitle(paste("Variable=", var_num, "Method=", cor_method, "Cell Num=", length(Cells(rna))))

dir.create(file.path(plot_dir, "umap"), showWarnings=FALSE)
ggsave(file.path(plot_dir, "umap",paste0(run.id, ".jpeg")))
# output coembed
coembed <- rbind(res$u, res$r)
dir.create(file.path(out_dir, "coembed"), showWarnings=FALSE)
write.csv(coembed, file=file.path(out_dir, "coembed", paste0(run.id, ".csv")), row.names = TRUE)

# transfer annotation with KNN classifier(k=30)
rna_cell_num <- ncol(rna)
rna_embed <- coembed[1:rna_cell_num, ]
atac_embed <- coembed[(rna_cell_num+1):nrow(coembed), ]
train_barcodes <- rownames(rna_embed)
train_annotations <- barcode_to_annotation_rna
train_annotations$barcodes <- paste0("RNA_", train_annotations$barcodes)
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
ggsave(file.path(plot_dir, "annotation_accuracy", paste0(run.id, ".jpeg")))
# save annotation result
cell_annotation_results <- data.frame(barcodes=rownames(atac_embed), predicted.id = ct.predictions, ct_truth=barcode_to_annotation_atac$ct_truth)
dir.create(file.path(out_dir, "ct_annotation"), showWarnings=FALSE)
write.csv(cell_annotation_results, file.path(out_dir, "ct_annotation", paste0(run.id, ".csv")))
   