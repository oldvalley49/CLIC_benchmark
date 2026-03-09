library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
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
# data_dir = directory for the RNA and ATAC data
# var_num = number of variable features to choose from
# total_num = number of total features to use for integration
# sub_num = subsampling size (in experiment, from 1 to 4 with 1 being 1/4 size and 4 being entire datasets)
# cor_method = method used to calculate conserved features; either bulk or pseudobulk
# index = index indicating the different runs under the same condition

   
# load data
rna_counts <- load_rna(data_dir)
atac_counts <- load_atac(data_dir)
if (activity_model == 'signac') {
    print('using signac gene activity')
    activity_counts = load_activity(data_dir)
} else {
    print('using maestro gene activity')
    activity_counts = load_maestro(data_dir)
}
barcode_to_annotation <- load_annotations(data_dir)

activity_counts <- activity_counts[intersect(rownames(activity_counts), 
                                             rownames(rna_counts)), ]
subsampled <- subsample(rna_counts, atac_counts, activity_counts, barcode_to_annotation, sub_num)
rna_counts <- subsampled$rna_counts
atac_counts <- subsampled$atac_counts
activity_counts <- subsampled$activity_counts
barcode_to_annotation <- subsampled$barcode_to_annotation

#initialize a Seurat object for RNA data and add celltype info
rna <- CreateSeuratObject(rna_counts, assay="RNA")
rna$ct_truth <- factor(barcode_to_annotation$ct_truth)
rna$orig.ident <- "RNA"
#initialize a Seurat object for ATAC data
granges.counts <- StringToGRanges(rownames(atac_counts))
granges.use <- seqnames(granges.counts) %in% standardChromosomes(granges.counts)
atac_counts <- atac_counts[as.vector(granges.use), ]

if (startsWith(tissue, "m")) {
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
} else{
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
}

seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- 'hg38'


chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    min.cells = 10,
    annotation = annotations
)
atac <- CreateSeuratObject(counts=chrom_assay, assay="ATAC")
atac$orig.ident <- "ATAC"
# initialize a Seurat object for ATAC activity data and add celltype info
atac[["ACTIVITY"]] <- CreateAssayObject(counts=activity_counts)
atac$ct_truth <- factor(barcode_to_annotation$ct_truth)

# standard preprocessing workflow for RNA
rna <- NormalizeData(rna)
if (startsWith(tissue, "m")) {
    out <- FindCLICFeatures(rna, score_name = paste0('mouse-', cor_method), initial_variable_features_num = var_num, nfeatures = total_num)
} else {
    out <- FindCLICFeatures(rna, score_name = paste0('human-', cor_method), initial_variable_features_num = var_num, nfeatures = total_num)
} 

use.features <- out$use_features
cor.num <- out$clic_num

# standard preprocessing workflow for ATAC
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)

# standard preprocessing for gene activity
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))

# generate run ID
run.id <-  paste("seurat", cor_method, var_num, cor.num, length(Cells(rna)), activity_model, index, sep="_")
print(run.id)
# identify anchors for cca
transfer.anchors <- FindTransferAnchors(reference=rna, query=atac, 
                                        features=use.features,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca")
# transfer cell type annotation
ct.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna$ct_truth,
                                    weight.reduction = atac[["lsi"]], dims = 2:30)
    
atac <- AddMetaData(atac, metadata = ct.predictions)
print(length(atac$ct_truth))
print(length(atac$predicted.id))
atac$annotation_correct <- atac$predicted.id == atac$ct_truth

# save annotation results 
cell_annotation_results <- data.frame(barcodes=Cells(atac), predicted.id = atac$predicted.id, ct_truth=atac$ct_truth)
dir.create(file.path("output/paired/rna_atac", tissue), showWarnings = FALSE)
dir.create(file.path("output/paired/rna_atac", tissue, "ct_annotation"), showWarnings = FALSE)
write.csv(cell_annotation_results, file.path("output/paired/rna_atac", tissue, "ct_annotation", paste0(run.id, ".csv")))
# plot annotation accuracy by cell type
predictions <- table(atac$ct_truth, atac$predicted.id)
predictions <- predictions/rowSums(predictions)  
predictions <- as.data.frame(predictions)

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + 
        geom_tile() + scale_fill_gradient(name = "Fraction of cells",
        low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") +
        ylab("Predicted cell type label (ATAC)") +
        theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle(paste("Variable=", var_num, "Method=", cor_method, "Cell Num=", length(Cells(rna))))

dir.create(file.path("plots/paired/rna_atac", tissue), showWarnings = FALSE)
dir.create(file.path("plots/paired/rna_atac", tissue, "annotation_accuracy"), showWarnings = FALSE)
ggsave(file.path("plots/paired/rna_atac", tissue, "annotation_accuracy", paste0(run.id, ".jpeg")) , width=20, height=15)

# coembed 
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
                            dims = 2:30)


atac[["RNA"]] <- imputation
coembed <- merge(x = rna, y = atac, add.cell.ids = c("RNA", "ATAC"))
coembed <- ScaleData(coembed, features = use.features, do.scale = FALSE)
coembed <- RunPCA(coembed, features = use.features, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

p2 <- DimPlot(coembed, group.by = c("orig.ident"))
dir.create(file.path("plots/paired/rna_atac", tissue, "umap"), showWarnings = FALSE)
ggsave(file.path("plots/paired/rna_atac", tissue, "umap",paste0(run.id, ".jpeg")))
dir.create(file.path("output/paired/rna_atac", tissue, "coembed"), showWarnings = FALSE)
write.csv(Embeddings(coembed, reduction = "pca"), 
            file = file.path("output/paired/rna_atac", tissue, "coembed", paste0(run.id, ".csv")), row.names = TRUE)

# coembed w/ cell type:
p3 <- DimPlot(coembed, group.by = "ct_truth", shuffle=TRUE)
dir.create(file.path("plots/paired/rna_atac", tissue, "umap_ct"), showWarnings = FALSE)
ggsave(file.path("plots/paired/rna_atac", tissue, "umap_ct",paste0(run.id, ".jpeg")))

# save use.features:
dir.create(file.path("output/paired/rna_atac", tissue, "features"), showWarnings = FALSE)
write.csv(use.features, file = file.path("output/paired/rna_atac", tissue, "features", paste0(run.id, ".csv")))


