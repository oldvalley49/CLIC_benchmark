library(Seurat)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(bindSC)
library(Matrix)
library(irlba)
library(umap)
library(Signac)
library(class)
library(cowplot)
library(GenomeInfoDb)
library(cluster)
library(glue)
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

# initialize seurat object RNA data and perform preprocessing including clustering step required for bindsc
rna <- CreateSeuratObject(rna_counts, assay="RNA")
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures=var_num)
rna <- ScaleData(rna)
rna <- RunPCA(rna)
rna <- FindNeighbors(rna, dims = 1:20, reduction="pca")
rna <- FindClusters(rna, resolution = 0.5)
colnames(rna) <- paste0("RNA_", colnames(rna))
# initialize seurat object for ATAC data and perform preprocessing including clustering step required for bindsc
colnames(atac_counts) <- paste0("ATAC_", colnames(atac_counts))
colnames(activity_counts) <- paste0("ATAC_", colnames(activity_counts))
granges.counts <- StringToGRanges(rownames(atac_counts))
granges.use <- seqnames(granges.counts) %in% standardChromosomes(granges.counts)
atac_counts <- atac_counts[as.vector(granges.use), ]
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
#seqlevelsStyle(annotations) <- 'UCSC'
#genome(annotations) <- "hg38"
chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c("-", "-"),
    genome = 'hg38',
    min.cells = 10,
    #annotation = annotations
)
atac <- CreateSeuratObject(counts=chrom_assay, assay="ATAC")
atac$orig.ident <- "ATAC"
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = "q0")
atac <- RunSVD(atac)
atac <- FindNeighbors(atac, dims = 1:20, reduction = "lsi")
atac <- FindClusters(atac, resolution = 0.5)
# filter genes
variable.features <- VariableFeatures(rna)
result <- feature_selection(variable.features, tissue, cor_method, var_num, total_num)
use.features <- result$use.features
cor.num <- result$cor.num

# generate run id
run.id <-  paste("bindsc", cor_method, var_num, cor.num, length(Cells(rna)), index, sep="_")
run.id

# we won't be able to use genes that don't have fragments in the promoter region and genomic region
use.features = intersect(rownames(activity_counts), use.features)

atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity_counts)
DefaultAssay(atac) <- "ACTIVITY"
atac <- NormalizeData(atac)

X <- rna[["RNA"]]$data[use.features, ]
Z0 <- atac[["ACTIVITY"]]$data[use.features,]
#Y <- atac_counts
# dimensional reduction
out <- dimReduce( dt1 =  X, dt2 = Z0,  K = 30)
x <- out$dt1
z0 <- out$dt2
y  <- atac@reductions$lsi@cell.embeddings
# run integration
res <- BiCCA( X = t(x) ,
                Y = t(y), 
                Z0 =t(z0), 
                X.clst = rna$seurat_clusters,
                Y.clst = atac$seurat_clusters,
                alpha = 0.5, 
                lambda = 0.5,
                K = 15,
                temp.path  = glue("bindsc-{tissue}-{var_num}-{cor_method}-{index}"),
                num.iteration = 50,
                tolerance = 0.01,
                save = TRUE,
                parameter.optimize = FALSE,
                block.size = 0)

# plot coembedded space
dim(res$u)
dim(res$r)
umap_plt <- umap(rbind(res$u, res$r))
umap_plt  <- data.frame("UMAP1"=umap_plt$layout[,1],
                        "UMAP2"=umap_plt$layout[,2],
                        "data" = c(rep("scRNA-seq",ncol(X)),
                                    rep("scATAC-seq",ncol(atac_counts))))

p <- ggplot(umap_plt, aes(x = UMAP1, y = UMAP2, color = data)) + 
    geom_point(size = 0.1) +
    theme_minimal() +
    ggtitle(paste("Variable=", var_num, "Method=", cor_method, "Cell Num=", length(Cells(rna))))
dir.create(file.path("plots/spatial", tissue), showWarnings=FALSE)
dir.create(file.path("plots/spatial", tissue, "umap"), showWarning=FALSE)
ggsave(file.path("plots/spatial", tissue, "umap",paste0(run.id, ".jpeg")), plot=p)
# output coembed
coembed <- rbind(res$u, res$r)
dir.create(file.path("output/spatial", tissue), showWarnings=FALSE)
dir.create(file.path("output/spatial", tissue, "coembed"), showWarnings=FALSE)
write.csv(coembed, file=file.path("output/spatial", tissue, "coembed", paste0(run.id, ".csv")), row.names = TRUE)

