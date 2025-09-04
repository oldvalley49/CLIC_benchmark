library(GenomicRanges)
library(glue)
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)
library(rtracklayer)

tissue <- "mBrain2"

rna_counts <- readMM("data/raw_data/rna_atac/spatial/mBrain2/5M_20um_RNA/matrix.mtx")
barcodes <- read.table(file="data/raw_data/rna_atac/spatial/mBrain2/5M_20um_RNA/barcodes.tsv", sep='\t')
genes <- read.table(file='data/raw_data/rna_atac/spatial/mBrain2/5M_20um_RNA/features.tsv', sep='\t')
colnames(rna_counts) <- barcodes$V1
rownames(rna_counts) <- genes$V2
rna_counts <- as(rna_counts, "CsparseMatrix")

frags <- "data/raw_data/rna_atac/spatial/mBrain2/GSM8494157_5M_20um_ATAC.fragments.tsv.gz"
fragments = CreateFragmentObject(frags)
peaks <- CallPeaks(fragments)
cells <- colnames(rna_counts)
atac_counts = FeatureMatrix(fragments = fragments, features=peaks)

rna_counts <- rna_counts[, colnames(atac_counts)]

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c("-", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


chrom_assay <- CreateChromatinAssay(
  counts=atac_counts,
  sep=c("-", "-"),
  fragments = frags,
  annotation = annotations
)

atac <- CreateSeuratObject(counts=chrom_assay)
gene.activities <- GeneActivity(atac, features = rownames(rna_counts))

rna_counts <- rowsum(as.matrix(rna_counts), group = rownames(rna_counts))
rna_counts <- as(rna_counts, "CsparseMatrix")
gene.activities <- rowsum(as.matrix(gene.activities), group=rownames(gene.activities))
gene.activities <- as(gene.activities, "CsparseMatrix")


out_dir <- "data/processed_data/spatial/mBrain2/"
writeMM(rna_counts, file = paste0(out_dir, "rna_counts.mtx"))
writeMM(gene.activities, file = paste0(out_dir, "activity_counts.mtx"))
writeMM(atac_counts, file = paste0(out_dir, "atac_counts.mtx"))
write.csv(colnames(rna_counts), file = paste0(out_dir, "barcode.csv"))
write.csv(rownames(rna_counts), file = paste0(out_dir, "genes.csv"))
write.csv(rownames(gene.activities), file=paste0(out_dir, "activity_genes.csv"))
write.csv(rownames(atac_counts), file=paste0(out_dir, "peaks.csv"))