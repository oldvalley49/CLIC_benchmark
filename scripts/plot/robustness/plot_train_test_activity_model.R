library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(glue)

# preprocess paired data
preprocess_paired <- function(results) {
  # validate input data
  required_cols <- c("var_num", "algorithm", "filter_method", "foscttm", "knn_auc", "asw", "ari", "asw.batch", "activity_model")
  missing_cols <- setdiff(required_cols, colnames(results))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # get baseline metrics
  baseline <- results[results$var_num == 2000, ]
  # get filtered metrics
  filtered <- results[results$var_num %in% c(4000, 5000,6000, 8000), ]
  # take average of metrics across replicates
  filtered <- filtered %>% 
    group_by(filter_method, activity_model) %>% 
    summarize(
      foscttm = mean(foscttm, na.rm = TRUE),
      knn_auc = mean(knn_auc, na.rm = TRUE),
      asw = mean(asw, na.rm = TRUE),
      ari = mean(ari, na.rm = TRUE),
      asw.batch = mean(asw.batch, na.rm = TRUE),
      .groups = "drop"
    )
              
   baseline <- baseline %>% 
    group_by(activity_model) %>% 
    summarize(
      foscttm = mean(foscttm, na.rm = TRUE),
      knn_auc = mean(knn_auc, na.rm = TRUE),
      asw = mean(asw, na.rm = TRUE),
      ari = mean(ari, na.rm = TRUE),
      asw.batch = mean(asw.batch, na.rm = TRUE),
      .groups = "drop"
    )
  
  baseline$filter_method <- "baseline"
  baseline$id <- glue("{baseline$filter_method}_{baseline$activity_model}")
  filtered$id <- glue("{filtered$filter_method}_{filtered$activity_model}")
  combined <- rbind(filtered, baseline)
  print(combined)
  return(combined)
}

# Process multiple datasets and average results
process_datasets <- function(file_list, preprocess_func) {
  # Process each dataset
  processed_datasets <- list()
  
  for (file_info in file_list) {
    tryCatch({
      data <- read.csv(file_info$path)
      processed_datasets[[file_info$name]] <- preprocess_func(data)
      message(paste("Successfully processed:", file_info$name))
    }, error = function(e) {
      warning(paste("Error processing", file_info$name, ":", e$message))
    })
  }
  print(processed_datasets)
  # Average the results
  
    avg_results <- bind_rows(processed_datasets, .id = "rep") %>%
        group_by(id) %>%
        summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
                    .groups = "drop")
  return(avg_results)
}

paired_datasets <- list(
  list(path = "output/paired/rna_atac/BMMC-s1d2/results.csv", name = "BMMC_s1d2"),
  list(path = "output/paired/rna_atac/BMMC-s4d8/results.csv", name = "BMMC_s4d8"),
  list(path = "output/paired/rna_atac/PBMC/results.csv", name = "PBMC"),
  list(path = "output/paired/rna_atac/TDBM/results.csv", name = "TDBM"),
  list(path = "output/paired/rna_atac/mBrain/results.csv", name = "mBrain"),
  list(path = "output/paired/rna_atac/mRetina/results.csv", name = "mRetina"),
  list(path = "output/paired/rna_atac/mSkin/results.csv", name = "mSkin")
)

paired_results <- process_datasets(paired_datasets, preprocess_paired)

maestro_results <- paired_results %>% filter(endsWith(id, "maestro"))
signac_results <- paired_results %>% filter(endsWith(id, "signac"))

signac_order <- c(
  "maestro-pearson_signac",
  "signac-pearson_signac",
  "baseline_signac"
)

maestro_order <- c(
  "maestro-pearson_maestro",
  "signac-pearson_maestro",
  "baseline_maestro"
)

signac_results <- signac_results %>%
  mutate(
    id = factor(id, levels = signac_order)
  )

maestro_results <- maestro_results %>%
  mutate(
    id = factor(id, levels = maestro_order)
  )
  
signac_results_scaled <- signac_results %>%
    mutate(across(
        where(is.numeric),
        ~ {
        rng <- range(.x, na.rm = TRUE)
        if (rng[1] == rng[2]) 0.5 else (.x - rng[1]) / (rng[2] - rng[1]) + 0.1
        }
    ))

maestro_results_scaled <- maestro_results %>%
    mutate(across(
        where(is.numeric),
        ~ {
        rng <- range(.x, na.rm = TRUE)
        if (rng[1] == rng[2]) 0.5 else (.x - rng[1]) / (rng[2] - rng[1]) + 0.1
        }
    ))

maestro_long <- maestro_results_scaled %>%
  pivot_longer(
    cols = -id,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("foscttm", "knn_auc", "asw", "ari", "asw.batch"))
  ) %>%
  mutate(id = factor(as.character(id), levels = maestro_order))


p <- ggplot(maestro_long, aes(x = metric, y = id)) +
  geom_point(aes(size = value, fill = value), shape = 21, color = "black") +
  scale_size(range = c(2, 14)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = "Paired results balloon plot",
    x = "Metric",
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

ggsave("plots/benchmark/train_test_activity/baloon_maestro.svg", width = 8, height = 3, p)

signac_long <- signac_results_scaled %>%
  pivot_longer(
    cols = -id,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("foscttm", "knn_auc", "asw", "ari", "asw.batch"))
  ) %>%
  mutate(id = factor(as.character(id), levels = signac_order))


p <- ggplot(signac_long, aes(x = metric, y = id)) +
  geom_point(aes(size = value, fill = value), shape = 21, color = "black") +
  scale_size(range = c(2, 14)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(
    title = "Paired results balloon plot",
    x = "Metric",
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

ggsave("plots/benchmark/train_test_activity/baloon_signac.svg", width = 8, height = 3, p)


