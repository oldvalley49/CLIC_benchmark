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

tissue <- "mEmbryo"

rna_fp <- "data/raw_data/rna_atac/spatial/processed/ME13_50um_50bc/GSM6799937_ME13_50um_matrix_merge.tsv.gz"
fragment_fp <- "data/raw_data/rna_atac/spatial/processed/ME13_50um_50bc/GSM6801813_ME13_50um_fragments.tsv.gz"

# get rna feature barcode matrix
rna_counts <- read.table(file = rna_fp, sep = '\t', header = TRUE)
rna_counts <- as.matrix(rna_counts, "dgCMatrix")
colnames(rna_counts) <- paste0(colnames(rna_counts), "-1")
rna_counts <- as(rna_counts, "CsparseMatrix")


peaks <- import("data/raw_data/rna_atac/spatial/processed/ME13_50um_50bc/CellRanger/GSM6801813/outs/peaks.bed")

# get spot barcodes
cells <- colnames(rna_counts)

# generate peak by spot matrix
frag = CreateFragmentObject(fragment_fp)
atac_counts = FeatureMatrix(fragments = frag, features = peaks, cells=cells)
atac_barcodes <- colnames(atac_counts)
all(cells==atac_barcodes)

# only take peaks that have more than 5 reads in total:
use.peaks <- rowSums(atac_counts) > 5
atac_counts <- atac_counts[use.peaks, ]

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c("-", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c("-", "-"),
  fragments = fragment_fp,
  annotation = annotations
)

atac <- CreateSeuratObject(counts=chrom_assay)
gene.activities <- GeneActivity(atac, features = rownames(rna_counts))

out_dir <- "data/processed_data/spatial/mEmbryo/"
writeMM(rna_counts, file = paste0(out_dir, "rna_counts.mtx"))
writeMM(gene.activities, file = paste0(out_dir, "activity_counts.mtx"))
writeMM(atac_counts, file = paste0(out_dir, "atac_counts.mtx"))
write.csv(cells, file = paste0(out_dir, "barcode.csv"))
write.csv(rownames(rna_counts), file = paste0(out_dir, "genes.csv"))
write.csv(rownames(gene.activities), file=paste0(out_dir, "activity_genes.csv"))
write.csv(rownames(atac_counts), file=paste0(out_dir, "peaks.csv"))
#write.csv(rna@meta.data[["ct3"]], file=paste0(out_dir, "ct_annotation.csv"))
















  