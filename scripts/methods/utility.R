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

load_annotations <- function(data_dir){
   # data frame containing correspondance from cell barcode to cell type
   ct_truth <- read.csv(file.path(data_dir, "ct_annotation.csv"))
   ct_truth <- ct_truth$x
   barcodes <- read.csv(file.path(data_dir, "barcode.csv"))
   barcodes <- barcodes$x
   barcode_to_annotation <- data.frame(
      barcodes = barcodes,
      ct_truth = ct_truth
   )
   return(barcode_to_annotation)
}

subsample <- function(rna_counts, atac_counts, activity_counts, barcode_to_annotation, sub_num){
   # filter cells with no activity

   cell_sums <- Matrix::colSums(activity_counts)
   nonempty_cells <- which(cell_sums > 0)
   nonempty_barcodes <- colnames(activity_counts)[nonempty_cells]
   # perform subsampling
   barcodes <- read.csv(file.path(data_dir, "barcode.csv"))
   barcodes <- barcodes$x
   barcodes <- intersect(barcodes, nonempty_barcodes)
   num_total <- length(barcodes)
   num_subsample <- floor(num_total*sub_num)
   barcode_subsampled <- sample(barcodes, num_subsample, replace=FALSE)
   rna_counts <- rna_counts[, barcode_subsampled]
   atac_counts <- atac_counts[, barcode_subsampled]
   activity_counts <- activity_counts[, barcode_subsampled]
   barcode_to_annotation <- barcode_to_annotation %>%
      dplyr::filter(barcodes %in% barcode_subsampled) %>%
      mutate(order = match(barcodes, barcode_subsampled)) %>%
      arrange(order) %>%
      dplyr::select(-order)
   print(dim(rna_counts))
   print(dim(atac_counts))
   print(dim(activity_counts))
   return(list(rna_counts=rna_counts, atac_counts=atac_counts, activity_counts=activity_counts, barcode_to_annotation=barcode_to_annotation))
}


subsample_no_annotation <- function(rna_counts, atac_counts, activity_counts, sub_num){
   # filter cells with no activity

   cell_sums <- Matrix::colSums(activity_counts)
   nonempty_cells <- which(cell_sums > 0)
   nonempty_barcodes <- colnames(activity_counts)[nonempty_cells]
   # perform subsampling
   barcodes <- read.csv(file.path(data_dir, "barcode.csv"))
   barcodes <- barcodes$x
   barcodes <- intersect(barcodes, nonempty_barcodes)
   num_total <- length(barcodes)
   num_subsample <- floor(num_total*sub_num)
   barcode_subsampled <- sample(barcodes, num_subsample, replace=FALSE)
   rna_counts <- rna_counts[, barcode_subsampled]
   atac_counts <- atac_counts[, barcode_subsampled]
   activity_counts <- activity_counts[, barcode_subsampled]
   print(dim(rna_counts))
   print(dim(atac_counts))
   print(dim(activity_counts))
   return(list(rna_counts=rna_counts, atac_counts=atac_counts, activity_counts=activity_counts))
}

filter_cor <- function(var_genes, cor_genes, n) {
   intersection <- c()  # Initialize an empty vector to store the intersection
   result <- c()        # Initialize an empty vector to store the result
   cor.num = 0
   for (gene_cor in cor_genes) {
      cor.num = cor.num + 1
      if (gene_cor %in% var_genes) {
         result <- c(result, gene_cor)
      }
      if (length(result) == n) {
         return(list(use.features = result, cor.num=cor.num))
         break
      }
   }
   
}

feature_selection <- function(variable.features, tissue, cor_method, var_num, total_num){
   if (startsWith(tissue, 'm')){
      spiecies <- "mouse"
   }
   else{
      spiecies <- "human"
   }
   corr <- read.csv(file.path("output/conserved_features", spiecies, paste0(cor_method, ".csv")))
   #order the genes by pearson correlation
   if (var_num!=total_num){
      corr$pearson_correlation <- as.numeric(corr$pearson_correlation)
      corr.sorted <- corr %>% arrange(desc(pearson_correlation))
      conserved.features <- corr.sorted$X
      use.features <- filter_cor(variable.features, conserved.features, total_num)
      cor.num <- use.features$cor.num
      use.features <- use.features$use.features
   } else{
      use.features <- variable.features
      cor.num <- Inf
   }
   return(list(use.features=use.features, cor.num=cor.num))
}

# given a coembedding, transfer labels using k-nearest neighbor classifier
transfer_label <- function(coembed, rna, barcode_to_annotation, k=30) {
   cell_num <- ncol(rna)
   rna_embed <- coembed[1:cell_num, ]
   atac_embed <- coembed[(cell_num+1):nrow(coembed), ]
   train_barcodes <- rownames(rna_embed)
   train_annotations <- barcode_to_annotation
   train_annotations$barcodes <- paste0("RNA_", train_annotations$barcodes)
   asserthtat(all(train_annotations$barcodes == train_barcodes), msg="cell barcodes are not concordant")
   ct.predictions <- knn(train = rna_embed, test = atac_embed, cl = train_annotations$ct_truth, k = k)
   return(list(ct.truth=train_annotations$ct_truth, ct.predictions = ct.predictions))
}

# peak.matrix is a data frame of counts (called peaks from atac raw data
# annotation file is a GFF3 file eg "Homo_sapiens.GRCh38.109.chr.gff3"
CreateGeneActivityMatrix <- function(
  peak.matrix,
  annotation.file,
  seq.levels = c(1:22, "X", "Y"),
  include.body = TRUE,
  upstream = 2000,
  downstream = 0,
  verbose = TRUE
) {
    # convert peak matrix to GRanges object
    peak.df <- rownames(x = peak.matrix)
    peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
    peak.df <- as.data.frame(x = peak.df)
    colnames(x = peak.df) <- c("chromosome", 'start', 'end')
    peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)

    # if any peaks start at 0, change to 1
    # otherwise GenomicRanges::distanceToNearest will not work
    BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
    # get annotation file, select genes
    #   anno <- rtracklayer::import(con = annotation.file)
    # anno <- GenomeInfoDb::keepSeqlevels(x = anno, value = seq.levels, pruning.mode = 'coarse')

    gtf <- rtracklayer::import(con = annotation.file)
    gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')

    # change seqlevelsStyle if not the same
    if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
        GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
    }
    gtf.genes <- gtf[gtf$type == 'gene']

    # 
    # Extend definition up/downstream
    if (include.body) {
        gtf.body_prom <- Signac::Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
    } else {
        gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
    }
    gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
    keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
    peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
    gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]

    # Some gtf rows will not have gene_name attribute
    # Replace it by gene_id attribute
    # 
    gene.ids$Name[is.na(gene.ids$Name)] <- gene.ids$gene_id[is.na(gene.ids$Name)]

    peak.ids$gene.name <- gene.ids$Name
    peak.ids <- as.data.frame(x = peak.ids)
    peak.ids$peak <- rownames(peak.matrix)[S4Vectors::queryHits(x = keep.overlaps)]
    annotations <- peak.ids[, c('peak', 'gene.name')]
    colnames(x = annotations) <- c('feature', 'new_feature')

    # collapse into expression matrix
    peak.matrix <- as(object = peak.matrix, Class = 'matrix')
    all.features <- unique(x = annotations$new_feature)

    mysapply <- if (verbose) pbapply::pbsapply else sapply
    
    newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x){

    features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature

    submat <- peak.matrix[features.use, ]

    if (length(x = features.use) > 1) {
        return(Matrix::colSums(x = submat))
    } else {
        return(submat)
    }
    })

    newmat <- t(x = newmat)

    rownames(x = newmat) <- all.features

    colnames(x = newmat) <- colnames(x = peak.matrix)

    return(as(object = newmat, Class = 'dgCMatrix'))
}