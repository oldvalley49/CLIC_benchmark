library(Seurat)
library(rliger)
library(dplyr)
library(class)
library(assertthat)
library(Matrix)

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