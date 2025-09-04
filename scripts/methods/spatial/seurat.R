library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(cowplot)
library(Matrix)
library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v86)
source("scripts/methods/utility.R")

# parse input argument
args <- commandArgs(trailingOnly = TRUE)
tissue = args[1]
var_num <- as.numeric(args[2])
total_num <- as.numeric(args[3])
sub_num <- as.numeric(args[4])
cor_method <- args[5]
index <- as.numeric(args[6])
data_dir <- paste0("data/processed_data/spatial/", tissue)
   
# load data
rna_counts <- load_rna(data_dir)
atac_counts <- load_atac(data_dir)
activity_counts <- load_activity(data_dir)

# subsample
subsampled <- subsample_no_annotation(rna_counts, atac_counts, activity_counts, sub_num)
rna_counts <- subsampled$rna_counts
atac_counts <- subsampled$atac_counts
activity_counts <- subsampled$activity_counts

#initialize a Seurat object for RNA data and add celltype info
rna <- CreateSeuratObject(rna_counts, assay="RNA")
rna$orig.ident <- "RNA"
#initialize a Seurat object for ATAC data
granges.counts <- StringToGRanges(rownames(atac_counts))
granges.use <- seqnames(granges.counts) %in% standardChromosomes(granges.counts)
atac_counts <- atac_counts[as.vector(granges.use), ]
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "hg38"
chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'hg38',
    min.cells = 10,
    # annotation = annotations
)
atac <- CreateSeuratObject(counts=chrom_assay, assay="ATAC")
atac$orig.ident <- "ATAC"
# initialize a Seurat object for ATAC activity data and add celltype info
atac[["ACTIVITY"]] <- CreateAssayObject(counts=activity_counts)

# standard preprocessing workflow for RNA
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures=var_num)
rna <- ScaleData(rna)
rna <- RunPCA(rna)

# standard preprocessing workflow for ATAC
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)

# standard preprocessing for gene activity
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)
atac <- ScaleData(atac, features = rownames(atac))

# filter the features
variable.features <- VariableFeatures(rna)
result <- feature_selection(variable.features, tissue, cor_method, var_num, total_num)
use.features <- result$use.features
cor.num <- result$cor.num
# generate run ID
run.id <-  paste("seurat", cor_method, var_num, cor.num, length(Cells(rna)), index, sep="_")
print(run.id)
# identify anchors for cca
transfer.anchors <- FindTransferAnchors(reference=rna, query=atac, 
                                        features=use.features,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca")

dir.create(file.path("output/spatial", tissue), showWarnings = FALSE)


dir.create(file.path("plots/spatial", tissue), showWarnings = FALSE)


# coembed 
refdata <- GetAssayData(rna, assay = "RNA", slot = "data")
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
                            dims = 2:30)
dir.create(file.path("output/spatial", tissue, "geximpute"))
gex.imputation.fp <- file.path("output/spatial", tissue, "geximpute", paste0(run.id, ".RDS"))
saveRDS(imputation@data, file=gex.imputation.fp)
atac[["RNA"]] <- imputation
coembed <- merge(x = rna, y = atac, add.cell.ids = c("RNA", "ATAC"))
coembed <- ScaleData(coembed, features = use.features, do.scale = FALSE)
coembed <- RunPCA(coembed, features = use.features, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

p2 <- DimPlot(coembed, group.by = c("orig.ident"))
dir.create(file.path("plots/spatial", tissue, "umap"), showWarnings = FALSE)
ggsave(file.path("plots/spatial", tissue, "umap",paste0(run.id, ".jpeg")))
dir.create(file.path("output/spatial", tissue, "coembed"), showWarnings = FALSE)
write.csv(Embeddings(coembed, reduction = "pca"), 
            file = file.path("output/spatial", tissue, "coembed", paste0(run.id, ".csv")), row.names = TRUE)

# impute ATAC profile

# transfer.anchors <- FindTransferAnchors(reference=atac, query=rna, 
#                                         features=use.features,
#                                         reference.assay = "ACTIVITY",
#                                         query.assay = "RNA",
#                                         reduction = "cca")

# refdata <- GetAssayData(atac, assay = "ATAC", slot = "data")

# atac_imputed <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = rna[["pca"]],
#                             dims = 1:30)
# # check if atac_imputed has rownames and column names
# if (startsWith(tissue, "m")){
#     gtf.fp <- "data/raw_data/gtf/mm10_2020-A.gtf"
# } else {
#     gtf.fp <- "data/raw_data/gtf/GRCh38_2020-A.gtf"
# }
# gene_level <- CreateGeneActivityMatrix(atac_imputed, annotation.file=gtf.fp)

# save as matrix file