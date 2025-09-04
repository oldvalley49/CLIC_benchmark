library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)
library(EnsDb.Hsapiens.v86)
source("scripts/methods/utility.R")




# runSeurat_unpaired <- function(rna_batch, atac_batch, rna_dir, atac_dir, var_num, total_num, sub_num, cor_method, index){

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
# # perform subsampling for RNA
# num_total <- length(barcodes_rna)
# num_subsample <- sub_num * num_total%/%4
# barcode_subsampled_rna <- sample(barcodes_rna, num_subsample, replace=FALSE)
# rna_counts <- rna_counts[, barcode_subsampled_rna]
# barcode_to_annotation_rna <- barcode_to_annotation_rna %>%
#    dplyr::filter(barcodes_rna %in% barcode_subsampled_rna) %>%
#    mutate(order = match(barcodes_rna, barcode_subsampled_rna)) %>%
#    arrange(order) %>%
#    dplyr::select(-order)
# # perform subsampling for ATAC
# num_total <- length(barcodes_atac)
# num_subsample <- sub_num * num_total%/%4
# barcode_subsampled_atac <- sample(barcodes_atac, num_subsample, replace=FALSE)
# atac_counts <- atac_counts[, barcode_subsampled_atac]
# activity_counts <- activity_counts[, barcode_subsampled_atac]
# barcode_to_annotation_atac <- barcode_to_annotation_atac %>%
#    dplyr::filter(barcodes_atac %in% barcode_subsampled_atac) %>%
#    mutate(order = match(barcodes_atac, barcode_subsampled_atac)) %>%
#    arrange(order) %>%
#    dplyr::select(-order)
#initialize a Seurat object for RNA data and add celltype info
rna <- CreateSeuratObject(rna_counts, assay="RNA")
rna$ct_truth <- factor(barcode_to_annotation_rna$ct_truth)
rna$orig.ident <- "RNA"
#initialize a Seurat object for ATAC data
# granges.counts <- StringToGRanges(rownames(atac_counts))
# granges.use <- seqnames(granges.counts) %in% standardChromosomes(granges.counts)
# atac_counts <- atac_counts[as.vector(granges.use), ]
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"
chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    # genome = 'hg38',
    # min.cells = 10,
    # annotation = annotations
)
atac <- CreateSeuratObject(counts=chrom_assay, assay="ATAC")
atac$orig.ident <- "ATAC"
# initialize a Seurat object for ATAC activity data and add celltype info
atac[["ACTIVITY"]] <- CreateAssayObject(counts=activity_counts)
atac$ct_truth <- factor(barcode_to_annotation_atac$ct_truth)

# standard preprocessing workflow for RNA
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures=var_num)
rna <- ScaleData(rna)
rna <- RunPCA(rna)

# standard preprocessing workflow for ATAC
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# standard preprocessing for gene activity
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))

# filter the features
variable.features <- VariableFeatures(rna)
result <- feature_selection(variable.features, rna_batch, cor_method, var_num, total_num)
use.features <- result$use.features
cor.num <- result$cor.num

# generate run ID
run.id <-  paste("seurat", rna_batch, atac_batch, cor_method, var_num, cor.num, length(Cells(rna)), length(Cells(atac)), index, sep="_")
print(run.id)
# identify anchors for ccaef
transfer.anchors <- FindTransferAnchors(reference=rna, query=atac, 
                                        features=use.features,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca")
# transfer cell type annotation
ct.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$ct_truth,
                                weight.reduction = atac[["lsi"]], dims = 2:30)
atac <- AddMetaData(atac, metadata = ct.predictions)
atac$annotation_correct <- atac$predicted.id == atac$ct_truth

# save annotation results
cell_annotation_results <- data.frame(barcodes=Cells(atac), predicted.id = atac$predicted.id, ct_truth=atac$ct_truth)
dir.create(file.path(out_dir, "ct_annotation"), showWarnings=FALSE)
write.csv(cell_annotation_results, file.path(out_dir, "ct_annotation", paste0(run.id, ".csv")))
# plot annotation accuracy by cell type
print(length(atac$ct_truth))
print(length(atac$predicted.id))
predictions <- table(atac$ct_truth, atac$predicted.id)
predictions <- predictions/rowSums(predictions)  
predictions <- as.data.frame(predictions)

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + 
    geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                    low = "#ffffc8", high = "#7d0025") + xlab("True cell type annotation (ATAC)") +
    ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle(paste("Variable=", var_num, "Method=", cor_method, "Cell Num=", length(Cells(atac))))
    
dir.create(file.path(plot_dir, "annotation_accuracy"), showWarnings=FALSE)
ggsave(file.path(plot_dir, "annotation_accuracy", paste0(run.id, ".jpeg")))

# coembed 
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")
refdata <- refdata[use.features, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
                            dims = 2:30)
dir.create(file.path("output/unpaired/rna_atac", paste0(rna_batch, "+", atac_batch), "geximpute"), showWarnings=FALSE)
gex.imputation.fp <- file.path("output/unpaired/rna_atac", paste0(rna_batch, "+", atac_batch), "geximpute", paste0(run.id, ".RDS"))
saveRDS(imputation@data, file=gex.imputation.fp)

atac[["RNA"]] <- imputation
coembed <- merge(x = rna, y = atac, add.cell.ids = c("RNA", "ATAC"))
coembed <- ScaleData(coembed, features = use.features, do.scale = FALSE)
coembed <- RunPCA(coembed, features = use.features, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

p2 <- DimPlot(coembed, group.by = c("orig.ident"))
dir.create(file.path(plot_dir, "umap"), showWarnings=FALSE)
ggsave(file.path(plot_dir, "umap",paste0(run.id, ".jpeg")))
dir.create(file.path(out_dir, "coembed"), showWarnings=FALSE)
write.csv(Embeddings(coembed, reduction = "pca"), 
            file = file.path(out_dir, "coembed", paste0(run.id, ".csv")), row.names = TRUE)

