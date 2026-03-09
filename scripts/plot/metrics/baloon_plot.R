library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(glue)


# min-max scale columns
min_max_scale <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# assign inverted_ranks (higher value means higher number)
inverse_rank <- function(x) {
  return(rank(x, ties.method = "min"))
}

# preprocess unpaired data
preprocess_unpaired <- function(results, filter_method) {

  # validate input data
  required_cols <- c("var_num", "algorithm", "filter_method", "asw", "ari", "asw.batch")
  missing_cols <- setdiff(required_cols, colnames(results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # get baseline metrics
  baseline <- results[results$var_num == 2000 & results$filter_method == "corr", ]
  # get filtered metrics
  filtered <- results[results$var_num %in% c(4000, 5000, 6000, 8000) & results$filter_method == "corr", ]
  filtered$filter_method <- "CLIC"
  # take average of metrics across replicates
  filtered <- filtered %>% 
    group_by(algorithm, filter_method, var_num) %>% 
    summarize(asw = mean(asw), ari = mean(ari), asw.batch = mean(asw.batch), .groups = "drop")
    
  baseline <- baseline %>% 
    group_by(algorithm) %>% 
    summarize(asw = mean(asw), ari = mean(ari), asw.batch=mean(asw.batch), .groups = "drop")

  baseline$filter_method <- "baseline"
  baseline$var_num <- 2000

  combined <- rbind(filtered, baseline)
  combined$id <- paste0(combined$filter_method, "_", combined$var_num)
  combined$filter_method <- NULL
  combined$var_num <- NULL
  
  # widen the data frame
  combined_wide <- combined %>%
    pivot_wider(
      names_from = algorithm, 
      values_from = c(asw, ari, asw.batch),
      names_glue = "{algorithm}_{.value}"
    )
  
  combined_wide <- as.data.frame(combined_wide)
  rownames(combined_wide) <- combined_wide$id
  combined_wide$id <- NULL
  
  return(combined_wide)
}

# preprocess paired data
preprocess_paired <- function(results, filter_method, activity_model) {
  # validate input data
  required_cols <- c("var_num", "algorithm", "filter_method", "foscttm", "knn_auc", "asw", "ari", "asw.batch", "activity_model")
  missing_cols <- setdiff(required_cols, colnames(results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # get baseline metrics
  baseline <- results[results$var_num == 2000 & results$filter_method == filter_method & results$activity_model == activity_model, ]
  # get filtered metrics
  filtered <- results[results$var_num %in% c(4000, 5000,6000, 8000) & results$filter_method == filter_method & results$activity_model == activity_model, ]
  filtered$filter_method <- "CLIC"
  # take average of metrics across replicates
  filtered <- filtered %>% 
    group_by(algorithm, filter_method, var_num) %>% 
    summarize(
      foscttm = mean(foscttm, na.rm = TRUE),
      knn_auc = mean(knn_auc, na.rm = TRUE),
      asw = mean(asw, na.rm = TRUE),
      ari = mean(ari, na.rm = TRUE),
      asw.batch = mean(asw.batch, na.rm = TRUE),
      .groups = "drop"
    )
              
  baseline <- baseline %>% 
    group_by(algorithm) %>% 
    summarize(
      foscttm = mean(foscttm, na.rm = TRUE),
      knn_auc = mean(knn_auc, na.rm = TRUE),
      asw = mean(asw, na.rm = TRUE),
      ari = mean(ari, na.rm = TRUE),
      asw.batch = mean(asw.batch, na.rm = TRUE),
      .groups = "drop"
    )
  
  baseline$filter_method <- "baseline"
  baseline$var_num <- 2000
  
  combined <- rbind(filtered, baseline)
  combined$id <- paste0(combined$filter_method, "_", combined$var_num)
  combined$filter_method <- NULL
  combined$var_num <- NULL
  
  # widen the data frame
  combined_wide <- combined %>%
    pivot_wider(
      names_from = algorithm, 
      values_from = c(foscttm, asw, ari, knn_auc, asw.batch),
      names_glue = "{algorithm}_{.value}"
    )
  
  combined_wide <- as.data.frame(combined_wide)
  rownames(combined_wide) <- combined_wide$id
  combined_wide$id <- NULL
  
  return(combined_wide)
}

# preprocess spatial data
preprocess_spatial <- function(results, filter_method) {
  # validate input data
  required_cols <- c("var_num", "algorithm", "filter_method", "foscttm", "knn_auc")
  missing_cols <- setdiff(required_cols, colnames(results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # get baseline metrics
  baseline <- results[results$var_num == 2000 & results$filter_method == "corr", ]
  # get filtered metrics
  filtered <- results[results$var_num %in% c(4000, 5000,6000, 8000) & results$filter_method == "corr", ]
  filtered$filter_method <- "CLIC"
  # take average of metrics across replicates
  filtered <- filtered %>% 
    group_by(algorithm, filter_method, var_num) %>% 
    summarize(foscttm = mean(foscttm), knn_auc = mean(knn_auc), .groups = "drop")
    
  baseline <- baseline %>% 
    group_by(algorithm) %>% 
    summarize(foscttm = mean(foscttm), knn_auc = mean(knn_auc), .groups = "drop")
  
  baseline$filter_method <- "baseline"
  baseline$var_num <- "2000"
  combined <- rbind(filtered, baseline)
  combined$id <- paste0(combined$filter_method, "_", combined$var_num)
  combined$filter_method <- NULL
  combined$var_num <- NULL
  
  # widen the data frame
  combined_wide <- combined %>%
    pivot_wider(
      names_from = algorithm, 
      values_from = c(foscttm, knn_auc),
      names_glue = "{algorithm}_{.value}"
    )
  
  combined_wide <- as.data.frame(combined_wide)
  rownames(combined_wide) <- combined_wide$id
  combined_wide$id <- NULL
  
  return(combined_wide)
}


# calculate rankings and normalize results
calculate_rankings <- function(avg_results) {
  ranking <- data.frame(avg_results)
  ranking <- ranking %>%
    mutate(across(where(is.numeric), inverse_rank, .names = "{.col}_inverse.rank"))
  
  ranking <- ranking %>%
    rowwise() %>%
    mutate(ranking = mean(c_across(ends_with("_inverse.rank")), na.rm = TRUE)) %>%
    ungroup()
    
  avg_results$average_inverted_ranking <- ranking$ranking
  
  # normalize the results
  avg_results <- avg_results %>%
    mutate(across(where(is.numeric), min_max_scale))
    
  return(avg_results)
}

# Create and save balloon plot
create_balloon_plot <- function(data, title, output_path, width, height) {
  data <- data[, colnames(data) != "average_inverted_ranking"]
  desired_order <- c("CLIC_8000", "CLIC_6000", "CLIC_5000", "CLIC_4000", "baseline_2000")
  data <- data[desired_order, ]
  p <- ggballoonplot(data, fill = "value") +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = title,
         x = "Metric",
         y = NA) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold", hjust=0.5)
    )
  
  # Create directory if it doesn't exist
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  
  # Save the plot
  ggsave(output_path, p, width = width, height = height)
}

# Process multiple datasets and average results
process_paired_datasets <- function(file_list, filter_method, activity_model, preprocess_func) {
  # Process each dataset
  processed_datasets <- list()
  
  for (file_info in file_list) {
    tryCatch({
      data <- read.csv(file_info$path)
      processed_datasets[[file_info$name]] <- preprocess_func(data, filter_method, activity_model)
      message(paste("Successfully processed:", file_info$name))
    }, error = function(e) {
      warning(paste("Error processing", file_info$name, ":", e$message))
    })
  }
  
  # Check if any datasets were successfully processed
  if (length(processed_datasets) == 0) {
    stop("No datasets were successfully processed.")
  }
  
  # Average the results
  avg_results <- Reduce("+", processed_datasets) / length(processed_datasets)
  
  return(avg_results)
}

# Process multiple datasets and average results
process_datasets <- function(file_list, preprocess_func) {
  # Process each dataset
  processed_datasets <- list()
  
  for (file_info in file_list) {
    tryCatch({
      data <- read.csv(file_info$path)
      processed_datasets[[file_info$name]] <- preprocess_func(data, filter_method)
      message(paste("Successfully processed:", file_info$name))
    }, error = function(e) {
      warning(paste("Error processing", file_info$name, ":", e$message))
    })
  }
  
  # Check if any datasets were successfully processed
  if (length(processed_datasets) == 0) {
    stop("No datasets were successfully processed.")
  }
  
  # Average the results
  avg_results <- Reduce("+", processed_datasets) / length(processed_datasets)
  
  return(avg_results)
}

create_bar_plot <- function(data, title, output_path) {
  desired_order <- c("baseline_2000", "CLIC_4000", "CLIC_5000", "CLIC_6000", "CLIC_8000")
  data$index <- factor(rownames(data), levels = desired_order)
  # print(data)
  p <- ggplot(data, aes(x = index, y = average_inverted_ranking)) +
    geom_bar(stat = "identity", fill = "steelblue", color='black', width = 0.6) +  # thinner bars
    theme_minimal() +
    labs(x = "Method", y = "Average Inverted Ranking", title = title) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),  # remove gridlines
    )
  ggsave(output_path, p, width = 4, height = 6)
}


# main analysis pipeline ----------------------------------------------------

# create base output directory if it doesn't exist
dir.create("plots/benchmark", showWarnings = FALSE, recursive = TRUE)

# Define dataset groups with file paths and names
unpaired_human_datasets <- list(
  list(path = "output/unpaired/rna_atac/BMMC-s1d2+BMMC-s2d4/results.csv", name = "s1d2_s2d4"), # similar size
  list(path = "output/unpaired/rna_atac/BMMC-s4d8+BMMC-s1d2/results.csv", name = "s4d8_s1d2"), # rna is larger
  list(path = "output/unpaired/rna_atac/BMMC-s3d6+BMMC-s2d4/results.csv", name = "s3d6_s2d4") # atac is larger
)

human_paired_datasets <- list(
  list(path = "output/paired/rna_atac/BMMC-s1d2/results.csv", name = "BMMC_s1d2"),
  list(path = "output/paired/rna_atac/BMMC-s4d8/results.csv", name = "BMMC_s4d8"),
  list(path = "output/paired/rna_atac/PBMC/results.csv", name = "PBMC"),
  list(path = "output/paired/rna_atac/TDBM/results.csv", name = "TDBM")
)

mouse_paired_datasets <- list(
  list(path = "output/paired/rna_atac/mBrain/results.csv", name = "mBrain"),
  list(path = "output/paired/rna_atac/mRetina/results.csv", name = "mRetina"),
  list(path = "output/paired/rna_atac/mSkin/results.csv", name = "mSkin")
)

human_spatial_datasets <- list(
  list(path = "output/spatial/brain/results.csv", name = "brain"),
  list(path = "output/spatial/Melanoma/results.csv", name = "Melanoma")
)


mouse_spatial_datasets <- list(
  list(path = "output/spatial/mEmbryo/results.csv", name = "mEmbryo"),
  list(path = "output/spatial/mEmbryo2/results.csv", name = "mEmbryo2"),
  list(path = "output/spatial/mBrain2/results.csv", name = "mBrain2")
)

# Process unpaired data
message("Processing unpaired data...")
unpaired_human_results <- process_datasets(unpaired_human_datasets, preprocess_unpaired)
unpaired_human_results <- calculate_rankings(unpaired_human_results)
create_balloon_plot(
  unpaired_human_results, 
  "BMMC Batch Effect Test", 
  "plots/benchmark/BMMCBatch_baloon.svg",
  8,
  4
)
create_bar_plot(
  unpaired_human_results,
  "BMMC Batch Effect",
  "plots/benchmark/BMMCBatch_bar.svg"
)

# Process paired human data
message("Processing paired human data...")
activity_models <- c('signac', 'maestro')
filter_methods <- c('maestro-pearson', 'signac-pearson')
for (activity_model in activity_models){
  for (filter_method in filter_methods){
    paired_human_results <- process_paired_datasets(human_paired_datasets, filter_method, activity_model, preprocess_paired)
    paired_human_results <- calculate_rankings(paired_human_results)
    create_balloon_plot(
      paired_human_results, 
      "Paired Human Datasets", 
      glue("plots/benchmark/paired_human_{filter_method}_{activity_model}_baloon.jpeg"),
      10,
      4
    )
    create_bar_plot(
      paired_human_results,
      "Paired Human Datasets",
      glue("plots/benchmark/paired_human_{filter_method}_{activity_model}_bar.jpeg")
    )
    create_balloon_plot(
      paired_human_results, 
      "Paired Human Datasets", 
      glue("plots/benchmark/paired_human_{filter_method}_{activity_model}_baloon.svg"),
      10,
      4
    )
    create_bar_plot(
      paired_human_results,
      "Paired Human Datasets",
      glue("plots/benchmark/paired_human_{filter_method}_{activity_model}_bar.svg")
    )
  }
}

# Process paired mouse data
message("Processing paired mouse data...")
activity_models <- c('maestro', 'signac')
filter_methods <- c('maestro-pearson', 'signac-pearson')

for (activity_model in activity_models){
  for (filter_method in filter_methods){
    paired_mouse_results <- process_paired_datasets(mouse_paired_datasets, filter_method, activity_model, preprocess_paired)
    paired_mouse_results <- calculate_rankings(paired_mouse_results)
    create_balloon_plot(
      paired_mouse_results, 
      "Paired Mouse Test", 
      glue("plots/benchmark/paired_mouse_{filter_method}_{activity_model}_baloon.jpeg"),
      10,
      4
    )
    create_bar_plot(
      paired_mouse_results,
      "Paired Mouse Datasets",
      glue("plots/benchmark/paired_mouse_{filter_method}_{activity_model}_bar.jpeg")
    )
    create_balloon_plot(
      paired_mouse_results, 
      "Paired Mouse Test", 
      glue("plots/benchmark/paired_mouse_{filter_method}_{activity_model}_baloon.svg"),
      10,
      4
    )
    create_bar_plot(
      paired_mouse_results,
      "Paired Mouse Datasets",
      glue("plots/benchmark/paired_mouse_{filter_method}_{activity_model}_bar.svg")
    )
  }
}



# Process spatial human data
message("Processing spatial human data...")
spatial_human_results <- process_datasets(human_spatial_datasets, preprocess_spatial)
spatial_human_results <- calculate_rankings(spatial_human_results)
create_balloon_plot(
  spatial_human_results, 
  "Spatial Human Test", 
  "plots/benchmark/spatial_human_baloon.svg",
  6,
  4
)
create_bar_plot(
  spatial_human_results,
  "Spatial Human Datasets",
  "plots/benchmark/spatial_human_bar.svg"
)

create_balloon_plot(
  spatial_human_results, 
  "Spatial Human Test", 
  "plots/benchmark/spatial_human_baloon.jpeg",
  6,
  4
)
create_bar_plot(
  spatial_human_results,
  "Spatial Human Datasets",
  "plots/benchmark/spatial_human_bar.jpeg"
)

# Process spatial mouse data
message("Processing spatial mouse data...")
spatial_mouse_results <- process_datasets(mouse_spatial_datasets, preprocess_spatial)
spatial_mouse_results <- calculate_rankings(spatial_mouse_results)
create_balloon_plot(
  spatial_mouse_results, 
  "Spatial Mouse Test", 
  "plots/benchmark/spatial_mouse_baloon.jpeg",
  6,
  4
)
create_bar_plot(
  spatial_mouse_results,
  "Spatial Mouse Datasets",
  "plots/benchmark/spatial_mouse_bar.jpeg"
)

create_balloon_plot(
  spatial_mouse_results, 
  "Spatial Mouse Test", 
  "plots/benchmark/spatial_mouse_baloon.jpeg",
  6,
  4
)
create_bar_plot(
  spatial_mouse_results,
  "Spatial Mouse Datasets",
  "plots/benchmark/spatial_mouse_bar.jpeg"
)


message("All analyses completed successfully!")