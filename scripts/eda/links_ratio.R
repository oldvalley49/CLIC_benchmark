library(Seurat)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(Signac)
library(data.table)
library(dplyr)

# testing overlap between SCT and normalize 
notSCT <- select(read.csv("data/processed_data/rna_atac/BMMC-s4d8/links_truth.csv"), gene, peak)
# SCT <- select(read.csv("data/processed_data/rna_atac/BMMC-s4d8/links_truth.csv"), gene, peak)
common <- inner_join(notSCT, SCT)

# notSCT = 17418 links
# SCT = 18091 links
# common = 13366 links


# taken from signac
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

# BMMC-s4d8 data


links <- select(read.csv("data/processed_data/rna_atac/BMMC-s4d8/links_truth.csv"), gene, peak)
links_gr <- makeGRangesFromDataFrame(links, keep.extra.columns=TRUE)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# unify transcripts
annotations_collapsed <- CollapseToLongestTranscript(annotations)

# get all possible pairs of peaks and genes that are 2000 bp close to each other
# see how many of peak-gene pairs are 2000 bp to each other -> how many of the assumptions were wrong?
# could also see how the including correlating features changes the results


truth <- select(read.csv("data/processed_data/rna_atac/BMMC-s4d8/links_truth.csv"), gene, peak)
truth_SCT <- dplyr::select(read.csv("data/processed_data/rna_atac/BMMC-s4d8/links_truth_SCT.csv"), gene, peak)
imputed <- dplyr::select(read.csv("/users/tfurutan/cross_modality/output/paired/rna_atac/BMMC-s4d8/imputed_links/seurat_pseudo_2000_Inf_9876_1.csv"), gene, peak)
imputed_8000 <- dplyr::select(read.csv("output/paired/rna_atac/BMMC-s4d8/imputed_links/seurat_pseudo_8000_5251_9876_1.csv"), gene, peak)
imputed_4000 <- dplyr::select(read.csv("output/paired/rna_atac/BMMC-s4d8/imputed_links/seurat_pseudo_4000_9280_9876_1.csv"), gene, peak)
common <- inner_join(truth, imputed)
common_SCT <- inner_join(truth_SCT, imputed)
common_8000 <- inner_join(truth, imputed_8000)
common_4000 <- inner_join(truth, imputed_4000)

# false positives
fp <- anti_join(imputed, truth)
fp_8000 <- anti_join(imputed_8000, truth)
fp_4000 <- anti_join(imputed_4000, truth)
