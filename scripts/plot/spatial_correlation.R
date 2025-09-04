library(ggplot2)
library(stringr)
library(glue)
library(dplyr)
library(Seurat)
source("scripts/methods/utility.R")

tissue <- "mEmbryo"
dir.path <- glue("output/spatial/{tissue}/impute_corr")
data_dir = "data/processed_data/spatial"
rna <- load_rna(file.path(data_dir, tissue))
rna <- CreateSeuratObject(rna, assay="RNA")
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, nfeatures=2000)
variable.features <- VariableFeatures(rna)

param_to_data <- list()

for (run.id in list.files(dir.path)) {
  data <- read.csv(glue("{dir.path}/{run.id}"), row.names = 1)
  data <- data[variable.features, ,drop = FALSE]
  run.id.parsed <- str_split(run.id, "_")[[1]]
  param <- run.id.parsed[3]
  if (!(as.numeric(param) >= 2000 && as.numeric(param) <= 8000)) {
    next
  }
  if (param %in% names(param_to_data)) {
    # merge existing and new data by row names
    param_to_data[[param]] <- merge(
      param_to_data[[param]], data, 
      by = "row.names", all = TRUE, suffixes = c("", "")
    )
    rownames(param_to_data[[param]]) <- param_to_data[[param]]$Row.names
    param_to_data[[param]]$Row.names <- NULL
  } else {
    param_to_data[[as.character(param)]] <- data
  }
}

summary_df <- data.frame(
  parameter = names(param_to_data),
  mean_cor_value = sapply(param_to_data, function(df) {
    mean(as.matrix(df), na.rm = TRUE)
  })
)

summary_df$parameter_numeric <- as.numeric(as.character(summary_df$parameter))
summary_df <- summary_df %>% arrange(parameter_numeric)
summary_df$parameter <- factor(summary_df$parameter, levels = summary_df$parameter)


p <- ggplot(summary_df, aes(x = parameter, y = mean_cor_value)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +  # thinner bars
  labs(title = "Average cor_value by Parameter",
       x = "Parameter", y = "Mean cor_value") +
  theme_minimal() +
  theme(
    axis.line = element_blank(),          # remove axis lines
    axis.ticks = element_blank(),         # remove axis ticks
    panel.grid = element_blank(),         # remove grid lines
    panel.background = element_blank(),   # remove panel background
    plot.background = element_blank()     # remove plot background
  )

ggsave(glue("plots/benchmark/spatial_{tissue}_imputation_correlation.svg"))


