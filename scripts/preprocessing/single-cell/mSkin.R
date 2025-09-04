library(dplyr)
library(Seurat)
library(Signac)
library(data.table)
library(tidyverse)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1234)

in_dir = "data/raw_data/rna_atac/mSkin"
celltype = read.csv(file.path(in_dir,"GSM4156597_skin_celltype.txt.gz"),sep='\t')

head(celltype)

frag_path = file.path(in_dir,"GSM4156597_skin.late.anagen.atac.fragments.bed.gz")
frg = data.table::fread(frag_path,sep = '\t',header = FALSE)
frg$V4 = stringr::str_replace_all(frg$V4, c(','='-','\\.'='-'))

cell_bc = stringr::str_replace_all(celltype$atac.bc, c(','='-','\\.'='-'))
rna_bc = stringr::str_replace_all(celltype$rna.bc, c(','='-','\\.'='-'))
fragments_bc_sel = frg %>% dplyr::filter(V4 %in% cell_bc)
fragments_bc_sel$V5 = rna_bc[match(fragments_bc_sel$V4,cell_bc)]

head(fragments_bc_sel)
dim(fragments_bc_sel)
head(fragments_bc_sel[,c(1:3,5)])
fragments_final = fragments_bc_sel[,c(1:3,5)]
fragments_final$V6=1
head(fragments_final)

in_dir = "data/raw_data/rna_atac/mSkin/"
cur_dir = getwd()
setwd(in_dir)

# intermediate output filename
out_path = "fragments_unsorted.tsv"
# final filename
filename = "mSkin_fragments.tsv"
pkg_source = "/jhpce/shared/libd/core/htslib/1.18/bin/"
data.table::fwrite(fragments_final,out_path,sep='\t',col.names=FALSE)

system(paste0("sort -k1,1 -k2,2n -k3,3n fragments_unsorted.tsv > ",filename),intern=TRUE)
system2(paste0(pkg_source, "bgzip"), paste0("-f ", filename))
system2(paste0(pkg_source, "tabix"), paste0("-p ", "bed ",filename,".gz"))
system("rm fragments_unsorted.tsv",intern=TRUE)

setwd(cur_dir)

rna_counts = fread(paste0(in_dir,"GSM4156608_skin.late.anagen.rna.counts.txt.gz"), sep="\t")

rna_genes <- rna_counts$gene
rna_counts <- rna_counts[,-1]
head(rna_counts)
rna_mat <- as(as.matrix(rna_counts), "sparseMatrix")
rownames(rna_mat)<-rna_genes
colnames(rna_mat)<-stringr::str_replace_all(colnames(rna_mat), c(','='-','\\.'='-'))

md<-read.csv(paste0(in_dir,"GSM4156597_skin_celltype.txt.gz"), sep="\t")
md <- md %>% 
    mutate(atac.bc = str_replace_all(atac.bc, c(','='-','\\.'='-')),
           rna.bc = str_replace_all(rna.bc, c(','='-','\\.'='-')),
           col = rna.bc) %>% 
    column_to_rownames("col")


atac_counts = Matrix::readMM(paste0(in_dir,"GSM4156597_skin.late.anagen.counts.txt.gz"))
bc = read.csv(paste0(in_dir,"GSM4156597_skin.late.anagen.barcodes.txt.gz"), sep="\t",header=FALSE)
head(bc)
peaks = read.csv(paste0(in_dir,"GSM4156597_skin.late.anagen.peaks.bed.gz"), sep="\t",header=FALSE)
peaks <- peaks %>% mutate(peaks = paste0(V1,"-",V2,"-",V3))
head(peaks)

colnames(atac_counts) <- str_replace_all(bc$V1, c(','='-','\\.'='-'))
rownames(atac_counts) <- peaks$peaks

cell_use = intersect(rownames(md),colnames(atac_counts))
length(cell_use)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
ucsc.levels <- stringr::str_replace(string=paste("chr",seqlevels(annotation),sep=""), pattern="chrMT", replacement="chrM")
seqlevels(annotation) <- ucsc.levels


atac_counts_fil = atac_counts[,cell_use]
atac_counts_fil = atac_counts_fil[rowSums(atac_counts_fil)>=1,]
chrom_assay <- CreateChromatinAssay(
    counts = atac_counts_fil,
    sep = c("-", "-"),
    fragments = file.path(in_dir, "mSkin_fragments.tsv.gz"),
    annotation = annotation
)
rna_mat = rna_mat[,cell_use]

seurat <- CreateSeuratObject(rna_mat[,cell_use],
                             assay = "RNA",
                             meta.data = md[cell_use,],
                             min.cells = 3)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^Mt-")
seurat[['ATAC']] <- chrom_assay

genes <- Features(seurat)
DefaultAssay(seurat) <- "ATAC"
gene.activities <- GeneActivity(seurat, features = genes)

out_dir <- "data/processed_data/rna_atac/mSkin/"
dir.create(out_dir, showWarnings = FALSE)
writeMM(seurat@assays[["RNA"]]@layers[["counts"]], file = paste0(out_dir, "rna_counts.mtx"))
writeMM(gene.activities, file = paste0(out_dir, "activity_counts.mtx"))
writeMM(seurat@assays[["ATAC"]]@counts, file = paste0(out_dir, "atac_counts.mtx"))
write.csv(Cells(seurat), file = paste0(out_dir, "barcode.csv"))
write.csv(genes, file = paste0(out_dir, "genes.csv"))
write.csv(rownames(gene.activities), file=paste0(out_dir, "activity_genes.csv"))
write.csv(Features(seurat), file=paste0(out_dir, "peaks.csv"))
write.csv(seurat@meta.data[["celltype"]], file=paste0(out_dir, "ct_annotation.csv"))



# ### for testing purposes, reconstrucinting seurat object and checking if clusters make sense
# # reconstruct seurat object
# data_dir = out_dir
# # read rna counts matrix as sparse matrix
# rna_counts_fp <- file.path(data_dir, "rna_counts.mtx")
# rna_counts <- readMM(rna_counts_fp)
# # add cell barcode and gene name to rna_counts
# barcodes <- read.csv(file.path(data_dir, "barcode.csv"))
# barcodes <- barcodes$x
# genes <- read.csv(file.path(data_dir, "genes.csv"))
# genes <- genes$x
# rownames(rna_counts) <- genes
# colnames(rna_counts) <- barcodes
# # read atac counts as sparse matrix
# atac_counts_fp <- file.path(data_dir, "atac_counts.mtx")
# atac_counts <- readMM(atac_counts_fp)
# peaks <- read.csv(file.path(data_dir, "peaks.csv"))
# peaks <- peaks$x
# rownames(atac_counts) <- peaks
# colnames(atac_counts) <- barcodes

# # read gene activity score 
# activity_counts_fp <- file.path(data_dir, "activity_counts.mtx")
# activity_counts <- readMM(activity_counts_fp)
# # add cell barcode and gene name to activity counts
# activity_genes <- read.csv(file.path(data_dir, "activity_genes.csv"))
# activity_genes <- activity_genes$x
# rownames(activity_counts) <- activity_genes
# colnames(activity_counts) <- barcodes

# # data frame containing correspondance from cell barcode to cell type
# ct_truth <- read.csv(file.path(data_dir, "ct_annotation.csv"))
# ct_truth <- ct_truth$x
# barcode_to_annotation <- data.frame(
#     barcodes = barcodes,
#     ct_truth = ct_truth
# )

# seurat <- CreateSeuratObject(rna_counts, assay="RNA")
# seurat$ct_truth <- factor(barcode_to_annotation$ct_truth)
# seurat$orig.ident <- "RNA"
# #initialize a Seurat object for ATAC data
# granges.counts <- StringToGRanges(rownames(atac_counts))
# granges.use <- seqnames(granges.counts) %in% standardChromosomes(granges.counts)
# atac_counts <- atac_counts[as.vector(granges.use), ]
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"
# chrom_assay <- CreateChromatinAssay(
#     counts = atac_counts,
#     sep = c(":", "-"),
#     genome = 'hg38',
#     min.cells = 10,
#     annotation = annotations
# )
# seurat[['ATAC']] <- chrom_assay
# seurat[["ACTIVITY"]] <- CreateAssayObject(counts=activity_counts)


# # check cell clustering
# DefaultAssay(seurat) <- "RNA"
# seurat <- SCTransform(seurat)
# seurat <- RunPCA(seurat)

# DefaultAssay(seurat) <- "ATAC"
# seurat <- FindTopFeatures(seurat, min.cutoff = 5)
# seurat <- RunTFIDF(seurat)
# seurat <- RunSVD(seurat)


# DefaultAssay(seurat) <- "ATAC"

# # first compute the GC content for each peak
# seurat <- RegionStats(seurat, genome = BSgenome.Hsapiens.UCSC.hg38)

# # link peaks to genes
# seurat <- LinkPeaks(
#    object = seurat,
#    peak.assay = "ATAC",
#    expression.assay = "SCT",
# )

# links <- Links(seurat)
# links.df <- as.data.frame(links)
# links.df <- links.df %>%
#    arrange(peak, -pvalue) %>%
#    filter(duplicated(peak) == FALSE) 

# write.csv(links.df, file = paste0(out_dir, "links_truth.csv"))
