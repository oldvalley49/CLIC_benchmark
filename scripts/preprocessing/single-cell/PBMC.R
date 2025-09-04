#future package max 
options(future.globals.maxSize = 8000 * 1024^2)

# Download libraries
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)

# Preprocess data according to the WNN tutorial from Seurat

in_dir = "data/raw_data/rna_atac/PBMC/"

inputdata.10x <- Read10X_h5(paste0(in_dir,"pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"))

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- paste0(in_dir,"pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz")
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts,
   sep = c(":", "-"),
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

pbmc <- subset(
   x = pbmc,
   subset = nCount_ATAC < 7e4 &
      nCount_ATAC > 5e3 &
      nCount_RNA < 25000 &
      nCount_RNA > 1000 &
      percent.mt < 20
)
# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# perform sub-clustering on cluster 6 to find additional structure
pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(pbmc) <- "sub.cluster"

# add annotations
pbmc <- RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
pbmc <- RenameIdents(pbmc, '0' = 'CD14 Mono', '9' ='CD14 Mono', '5' = 'CD16 Mono')
pbmc <- RenameIdents(pbmc, '10' = 'Naive B', '11' = 'Intermediate B', '17' = 'Memory B', '21' = 'Plasma')
pbmc <- RenameIdents(pbmc, '7' = 'NK')
pbmc <- RenameIdents(pbmc, '4' = 'CD4 TCM', '13'= "CD4 TEM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
pbmc <- RenameIdents(pbmc, '2' = 'CD8 Naive', '8'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2', '6_4' ='CD8 TEM_2')
pbmc <- RenameIdents(pbmc, '18' = 'MAIT')
pbmc <- RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')
pbmc$ct1 <- Idents(pbmc)


pbmc@meta.data$ct2 <- pbmc$ct1

#merge annotations by similar cell type
levels(pbmc@meta.data$ct2) <- c("other_T","other_T","CD8 T","CD8 T","CD8 T","CD4 T","CD4 T","CD4 T",
                                "CD4 T","NK","B","B","B","other","Mono","Mono","DC","other","DC")
pbmc@meta.data$ct3 <- pbmc$ct1

levels(pbmc@meta.data$ct3) <- c("other_T","other_T","CD8 Naive","CD8 TEM","CD8 TEM","CD4 T","CD4 T","CD4 T",
                                "CD4 T","NK","B","B","B","other","Mono","Mono","DC","other","DC")
#w e use ct3 for annotation
pbmc$ct1 = as.character(pbmc$ct1)
pbmc$ct2 = as.character(pbmc$ct2)
pbmc$ct3 = as.character(pbmc$ct3)


p1 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "ct1", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p2 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "ct2", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "ct3", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

combined_plot <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

# save plot to file
ggsave("plots/preprocess/10x_pbmc_annotation.png", plot = combined_plot, width = 10, height = 4, dpi = 300)

# pbmc$assign_a <- sample(x=1:4, size=ncol(pbmc), replace=TRUE)
# pbmc$assign_b <- sample(x=1:4, size=ncol(pbmc), replace=TRUE)
# pbmc$assign_c <- sample(x=1:4, size=ncol(pbmc), replace=TRUE)

# split object into rna and atac
DefaultAssay(pbmc) <- "RNA"
rna <- DietSeurat(pbmc, assays="RNA")
DefaultAssay(pbmc) <- "ATAC"
atac <- DietSeurat(pbmc, assays="ATAC")

# quantify gene activity
gene.activities <- GeneActivity(atac, features = Features(rna))

# export dat
out_dir <- "data/processed_data/rna_atac/PBMC/"
writeMM(rna@assays[["RNA"]]@layers[["counts"]], file = paste0(out_dir, "rna_counts.mtx"))
writeMM(gene.activities, file = paste0(out_dir, "activity_counts.mtx"))
writeMM(atac@assays[["ATAC"]]@counts, file = paste0(out_dir, "atac_counts.mtx"))
write.csv(Cells(rna), file = paste0(out_dir, "barcode.csv"))
write.csv(Features(rna), file = paste0(out_dir, "genes.csv"))
write.csv(rownames(gene.activities), file=paste0(out_dir, "activity_genes.csv"))
write.csv(Features(atac), file=paste0(out_dir, "peaks.csv"))
write.csv(rna@meta.data[["ct3"]], file=paste0(out_dir, "ct_annotation.csv"))

#a$x


