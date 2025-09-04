library(Matrix)
library(stringr)
library(glue)



args <- commandArgs(trailingOnly=TRUE)
tissue = args[1]
method = args[2]
method_name = str_split(method, "\\.")
method_name = method_name[[1]][1]

truth <- readRDS(glue("data/processed_data/rna_atac/{tissue}/norm_rna_counts.RDS"))
imputed <- readRDS(glue("output/paired/rna_atac/{tissue}/geximpute/{method}"))

rowwise_pearson <- function(x, y) {
  x_centered <- x - rowMeans(x)
  y_centered <- y - rowMeans(y)

  numerator <- rowSums(x_centered * y_centered)
  denominator <- sqrt(rowSums(x_centered^2)) * sqrt(rowSums(y_centered^2))

  numerator / denominator
}

genes <- rownames(truth)
cells <- colnames(truth)

imputed <- imputed[genes, cells, drop=FALSE]
correlation <- rowwise_pearson(truth, imputed)

correlation_df <- data.frame(
    row.names = genes,
    cor_value = correlation
)

write.csv(correlation_df, file=glue("output/spatial/{tissue}/impute_corr/{method_name}.csv"))

