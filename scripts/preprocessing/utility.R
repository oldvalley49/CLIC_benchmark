library(Seurat)
library(rliger)
library(dplyr)
library(class)
library(assertthat)
library(Matrix)
library(glue)

load_rna <- function(data_dir) {
   rna_counts_fp <- file.path(data_dir, "rna_counts.mtx")
   rna_counts <- readMM(rna_counts_fp)
   # add cell barcode and gene name to rna_counts
   barcodes <- read.csv(file.path(data_dir, "barcode.csv"))
   barcodes <- barcodes$x
   genes <- read.csv(file.path(data_dir, "genes.csv"))
   genes <- genes$x
   rownames(rna_counts) <- genes
   colnames(rna_counts) <- barcodes
   return(rna_counts)
}

load_atac <- function(data_dir){
   # read atac counts as sparse matrix
   atac_counts_fp <- file.path(data_dir, "atac_counts.mtx")
   atac_counts <- readMM(atac_counts_fp)
   barcodes <- read.csv(file.path(data_dir, "barcode.csv"))
   barcodes <- barcodes$x
   peaks <- read.csv(file.path(data_dir, "peaks.csv"))
   peaks <- peaks$x
   rownames(atac_counts) <- peaks
   colnames(atac_counts) <- barcodes
   return(atac_counts)
}

load_activity <- function(data_dir){
   # read gene activity score 
   activity_counts_fp <- file.path(data_dir, "activity_counts.mtx")
   activity_counts <- readMM(activity_counts_fp)
   # add cell barcode and gene name to activity counts
   activity_genes <- read.csv(file.path(data_dir, "activity_genes.csv"))
   activity_genes <- activity_genes$x
   barcodes <- read.csv(file.path(data_dir, "barcode.csv"))
   barcodes <- barcodes$x
   rownames(activity_counts) <- activity_genes
   colnames(activity_counts) <- barcodes
   return(activity_counts)
}

output_BED <- function(data_dir, peaks){
   coords <- do.call(rbind, strsplit(peaks, "-"))
   bed <- data.frame(
      chrom = coords[,1],
      start = as.integer(coords[,2]),
      end   = as.integer(coords[,3])
   )
   write.table(
      bed,
      file = file.path(data_dir, "peaks.bed"),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
   )
}

output_TSV <- function(data_dir, barcodes){
   write.table(
      barcodes,
      file = file.path(data_dir, "barcodes.tsv"),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
   )
}

compute_MAESTRO <- function(data_dir, genome){ # genome: GRCh38 or GRCm38

   barcodes <- read.csv(file.path(data_dir, "barcode.csv"))
   barcodes <- barcodes$x
   peaks <- read.csv(file.path(data_dir, "peaks.csv"))
   peaks <- peaks$x

   output_BED(data_dir, peaks)
   output_TSV(data_dir, barcodes)

   # MAESTRO
   cmd <- glue(
   "conda run -n MAESTRO MAESTRO scatac-genescore ",
   "--format mtx ",
   "--peakcount {data_dir}/atac_counts.mtx ",
   "--feature {data_dir}/peaks.bed ",
   "--barcode {data_dir}/barcodes.tsv ",
   "--species {genome} ",
   "--model Enhanced ",
   "-d {data_dir} ",
   "--outprefix MAESTRO"
   )

   cat(cmd, "\n")
   system(cmd)

   rna_counts <- load_rna(data_dir)
   gene.activities <- Read10X_h5(glue('{data_dir}/MAESTRO_gene_score.h5'))
   gene.activities <- gene.activities[intersect(rownames(rna_counts), rownames(gene.activities)), colnames(rna_counts)]
   writeMM(gene.activities, file = glue('{data_dir}/maestro_counts.mtx'))
   write.csv(rownames(gene.activities), file=glue("{data_dir}/maestro_genes.csv"))
}