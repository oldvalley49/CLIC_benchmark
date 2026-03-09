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
set.seed(1234)

in_dir = "data/raw_data/rna_atac/mBrain"
celltype = read.csv(file.path(in_dir,"GSM4156599_brain_celltype.txt.gz"),sep='\t')

head(celltype)

frag_path = file.path(in_dir,"GSM4156599_brain.atac.fragments.bed.gz")
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

in_dir = "data/raw_data/rna_atac/mBrain/"
cur_dir = getwd()
setwd(in_dir)

# intermediate output filename
out_path = "fragments_unsorted.tsv"
# final filename
filename = "mBrain_fragments.tsv"
pkg_source = "/jhpce/shared/libd/core/htslib/1.18/bin/"
data.table::fwrite(fragments_final,out_path,sep='\t',col.names=FALSE)

system(paste0("sort -k1,1 -k2,2n -k3,3n fragments_unsorted.tsv > ",filename),intern=TRUE)
system2(paste0(pkg_source, "bgzip"), paste0("-f ", filename))
system2(paste0(pkg_source, "tabix"), paste0("-p ", "bed ",filename,".gz"))
system("rm fragments_unsorted.tsv",intern=TRUE)

setwd(cur_dir)

rna_counts = fread(paste0(in_dir,"GSM4156610_brain.rna.counts.txt.gz"), sep="\t")

rna_genes <- rna_counts$gene
rna_counts <- rna_counts[,-1]
head(rna_counts)
rna_mat <- as(as.matrix(rna_counts), "sparseMatrix")
rownames(rna_mat)<-rna_genes
colnames(rna_mat)<-stringr::str_replace_all(colnames(rna_mat), c(','='-','\\.'='-'))

md<-read.csv(paste0(in_dir,"GSM4156599_brain_celltype.txt.gz"), sep="\t")
md <- md %>% 
    mutate(atac.bc = str_replace_all(atac.bc, c(','='-','\\.'='-')),
           rna.bc = str_replace_all(rna.bc, c(','='-','\\.'='-')),
           col = rna.bc) %>% 
    column_to_rownames("col")


atac_counts = Matrix::readMM(paste0(in_dir,"GSM4156599_brain.counts.txt.gz"))
bc = read.csv(paste0(in_dir,"GSM4156599_brain.barcodes.txt.gz"), sep="\t",header=FALSE)
head(bc)
peaks = read.csv(paste0(in_dir,"GSM4156599_brain.peaks.bed.gz"), sep="\t",header=FALSE)
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
    fragments = file.path(in_dir, "mBrain_fragments.tsv.gz"),
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

out_dir <- "data/processed_data/rna_atac/mBrain/"
dir.create(out_dir, showWarnings = FALSE)
writeMM(seurat@assays[["RNA"]]@layers[["counts"]], file = paste0(out_dir, "rna_counts.mtx"))
writeMM(gene.activities, file = paste0(out_dir, "activity_counts.mtx"))
writeMM(seurat@assays[["ATAC"]]@counts, file = paste0(out_dir, "atac_counts.mtx"))
write.csv(Cells(seurat), file = paste0(out_dir, "barcode.csv"))
write.csv(genes, file = paste0(out_dir, "genes.csv"))
write.csv(rownames(gene.activities), file=paste0(out_dir, "activity_genes.csv"))
write.csv(Features(seurat), file=paste0(out_dir, "peaks.csv"))
write.csv(seurat@meta.data[["celltype"]], file=paste0(out_dir, "ct_annotation.csv"))



