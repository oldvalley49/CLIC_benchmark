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
library(IRanges)
library(data.table)

# processing
matrix <- readMM("data/raw_data/rna_atac/spatial/Melanoma/expression/6426588be2c9c436276e802d/matrix.mtx")
barcode <- read.table("data/raw_data/rna_atac/spatial/Melanoma/expression/6426588be2c9c436276e802d/barcodes.tsv", sep='\t', header=FALSE)
features <- read.table("data/raw_data/rna_atac/spatial/Melanoma/expression/6426588be2c9c436276e802d/features.tsv", sep='\t', header=FALSE)
rownames(matrix) <- features$V2
colnames(matrix) <- barcode$V1

# split into RNA and ATAC
rna <- matrix[which(features$V3=='Gene Expression'), ]
atac <- matrix[which(features$V3=='Peaks'), ]
rownames(atac) <- gsub(":", "-", rownames(atac))


# generate activity from peak by cell

peak_coords <- do.call(rbind, strsplit(rownames(atac), "[-]"))
peak_gr <- GRanges(seqnames=peak_coords[,1],
                    ranges = IRanges(start = as.integer(peak_coords[, 2]),
                    end = as.integer(peak_coords[,3]))
                    )

# load annotation
annotation <- import("data/raw_data/gtf/GRCh38_2020-A.gtf")
annotation <- annotation[annotation$type=='gene', ]
annotation <- Extend(annotation, upstream=2000)

annotation <- annotation[annotation$gene_name %in% rownames(rna)]
annotation <- annotation[match(rownames(rna), annotation$gene_name)]


# find overlaps
hits <- findOverlaps(peak_gr, annotation)
peak_to_gene <- sparseMatrix(
  i = subjectHits(hits),
  j = queryHits(hits),
  x = 1,
  dims = c(length(annotation), length(peak_gr)),
  dimnames = list(annotation$gene_name, rownames(atac))
)
activity <- peak_to_gene %*% atac


rna <- rowsum(as.matrix(rna), group = rownames(rna))
rna <- as(rna, "CsparseMatrix")
activity <- rowsum(as.matrix(activity), group=rownames(activity))
activity <- as(activity, "CsparseMatrix")

out_dir <- "data/processed_data/spatial/Melanoma/"
writeMM(rna, file = paste0(out_dir, "rna_counts.mtx"))
writeMM(activity, file = paste0(out_dir, "activity_counts.mtx"))
writeMM(atac, file = paste0(out_dir, "atac_counts.mtx"))
write.csv(colnames(rna), file = paste0(out_dir, "barcode.csv"))
write.csv(rownames(rna), file = paste0(out_dir, "genes.csv"))
write.csv(rownames(activity), file=paste0(out_dir, "activity_genes.csv"))
write.csv(rownames(atac), file=paste0(out_dir, "peaks.csv"))