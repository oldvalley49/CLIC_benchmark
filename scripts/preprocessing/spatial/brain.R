library(GenomicRanges)
library(glue)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(Matrix)
library(rtracklayer)
# make sure to activate environment for macs2

tissue <- "brain"

# for some reason could'nt find the matrix file for rna in wanlu's directory so did it locally on my laptop w/ the copy i got beforehand
#rna_fp <- "data/raw_data/rna_atac/spatial/processed/HumanBrain_50um_50bc/GSM6206885_HumanBrain_50um_matrix.tsv"

fragment_fp <- "data/raw_data/rna_atac/spatial/processed/HumanBrain_50um_50bc/GSM6206884_HumanBrain_50um_fragments.tsv.gz"
rna_counts <- readMM("data/processed_data/spatial/brain/rna_counts.mtx")
barcode <- read.csv("data/processed_data/spatial/brain/barcode.csv")
barcode <- barcode$x
genes <- read.csv("data/processed_data/spatial/brain/genes.csv")
genes <- genes$x

colnames(rna_counts) <- barcode
rownames(rna_counts) <- genes

# get peaks
peaks <- import("data/raw_data/rna_atac/spatial/processed/HumanBrain_50um_50bc/peaks.bed")

# get reference peaks and filter blacklist regions
# load("data/processed_data/references/ENCODE_human_mouse_DHS_noY.rda")
# rm(mouse_regions)
# load("data/processed_data/references/hg38-blacklist.v2.rda")
# ref_peaks_filtered = human_regions[human_regions %outside% human_blacklist_granges,]
# # merge overlapping regions
# ref_peaks_filtered = GenomicRanges::reduce(ref_peaks_filtered)

# get spot barcodes
cells <- colnames(rna_counts)

# generate peak by spot matrix 
frag = CreateFragmentObject("data/raw_data/rna_atac/spatial/processed/HumanBrain_50um_50bc/GSM6206884_HumanBrain_50um_fragments.tsv.gz")
# ref_peaks_filtered_dummy <- data.frame((ref_peaks_filtered_dummy))
atac_counts = FeatureMatrix(fragments = frag, features = peaks, cells=cells)

# create count matrix
# atac_counts <- rsparsematrix(5, length(cells), nnz=5)
# rownames(atac_counts) <- paste0(ref_peaks_filtered_dummy$seqnames, ":", ref_peaks_filtered_dummy$start, "-", ref_peaks_filtered_dummy$end)
colnames(atac_counts) <- cells
atac_barcodes <- colnames(atac_counts)


# only take peaks that have more than 5 reads in total:
use.peaks <- rowSums(atac_counts) > 5
atac_counts <- atac_counts[use.peaks, ]

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c("-", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragment_fp,
  annotation = annotations
)

atac <- CreateSeuratObject(counts=chrom_assay)

gene.activities <- GeneActivity(atac, features = rownames(rna_counts))

out_dir <- "data/processed_data/spatial/brain/"
writeMM(rna_counts, file = paste0(out_dir, "rna_counts.mtx"))
writeMM(gene.activities, file = paste0(out_dir, "activity_counts.mtx"))
writeMM(atac_counts, file = paste0(out_dir, "atac_counts.mtx"))
write.csv(cells, file = paste0(out_dir, "barcode.csv"))
write.csv(rownames(rna_counts), file = paste0(out_dir, "genes.csv"))
write.csv(rownames(gene.activities), file=paste0(out_dir, "activity_genes.csv"))
write.csv(rownames(atac_counts), file=paste0(out_dir, "peaks.csv"))

















  