library(tidyverse)

# Input / Output Setup
out_dir <- "addendum/res"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ONT files
file_paths_ONT <- list.files(path = "data/nanopore", pattern = "nano_.*\\.rds$", full.names = TRUE)
file_order_ONT <- order(as.numeric(str_extract(basename(file_paths_ONT), "\\d+")))
file_paths_ONT <- file_paths_ONT[file_order_ONT]

# EM files
file_paths_EM <- list.files(path = "data/em-seq", pattern = "em_.*\\.rds$", full.names = TRUE)
file_order_EM <- order(as.numeric(str_extract(basename(file_paths_EM), "\\d+")))
file_paths_EM <- file_paths_EM[file_order_EM]

# Initialize empty lists to store results
results_list <- list()
conf_mat_list <- list() # Added list for confusion matrices

# Helper function to calculate Positive Agreement (Intersection over Union)
calc_positive_agreement <- function(mat, category) {
  intersection <- mat[category, category]
  union <- sum(mat[category, ]) + sum(mat[, category]) - intersection
  if (union == 0) return(NA)
  return(intersection / union)
}

# OUTER LOOP: Iterate over samples
for (i in seq_along(file_paths_ONT)) {
  
  sample_id <- str_extract(basename(file_paths_ONT[i]), "\\d+")
  message(sprintf("Loading and processing Sample: %s...", sample_id))
  
  # Load matched data into memory
  ont_data <- readRDS(file_paths_ONT[i])
  em_data  <- readRDS(file_paths_EM[i])
  
  # Adjust positions to match
  ont_data <- ont_data |>
    mutate(pos = paste(chr, start+1, sep = "_")) |> 
    filter(score >= 10)
  
  em_data  <- em_data  |> 
    mutate(pos = paste(chr, start, sep = "_")) |> 
    filter(score >= 10)
  
  # Join directly into a wide format and sample 1M shared loci
  shared_wide <- inner_join(
    ont_data %>% select(pos, beta_ONT = beta), 
    em_data %>% select(pos, beta_EM = beta), 
    by = "pos"
  ) 
  
  # Bin both technologies independently
  binned_data <- shared_wide %>%
    mutate(
      bin_EM = factor(cut(beta_EM, breaks = c(-Inf, 0.2, 0.8, Inf), labels = c("Low", "Medium", "High")), 
                      levels = c("Low", "Medium", "High")),
      bin_ONT = factor(cut(beta_ONT, breaks = c(-Inf, 0.2, 0.8, Inf), labels = c("Low", "Medium", "High")), 
                       levels = c("Low", "Medium", "High"))
    )
  
  # Create the 3x3 confusion matrix (Added names for clear axes)
  conf_mat <- table(EM = binned_data$bin_EM, ONT = binned_data$bin_ONT)
  
  # Calculate overall accuracy and category-specific positive agreement
  overall_acc <- sum(diag(conf_mat)) / sum(conf_mat)
  pa_low      <- calc_positive_agreement(conf_mat, "Low")
  pa_med      <- calc_positive_agreement(conf_mat, "Medium")
  pa_high     <- calc_positive_agreement(conf_mat, "High")
  
  # Print quick summary to console
  message(sprintf("  Overall Acc: %.1f%% | PA Low: %.1f%% | PA Med: %.1f%% | PA High: %.1f%%", 
                  overall_acc * 100, pa_low * 100, pa_med * 100, pa_high * 100))
  
  # Store results in the lists
  results_list[[i]] <- tibble(
    Sample_ID = sample_id,
    Overall_Accuracy = overall_acc,
    PA_Low = pa_low,
    PA_Medium = pa_med,
    PA_High = pa_high
  )
  
  # Store the confusion matrix using the sample_id as the name
  conf_mat_list[[sample_id]] <- conf_mat
  
  # UNLOAD: Free up memory before loading the next sample
  rm(ont_data, em_data, shared_wide, binned_data, conf_mat)
  gc() 
}

# Bind all results together and save as a CSV
final_results_df <- bind_rows(results_list) |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(final_results_df, file.path(out_dir, "positive_agreement_summary_EM.csv"))


# --- Process and Save Confusion Matrices ---

# 1. Flatten the list of matrices into a single dataframe
formatted_matrices <- bind_rows(lapply(names(conf_mat_list), function(sid) {
  mat <- conf_mat_list[[sid]]
  # Convert table to dataframe
  df <- as.data.frame(as.table(mat)) 
  df$Sample <- sid
  return(df)
}))

# 2. Reshape into a publication-friendly wide format
supp_table_df <- formatted_matrices |>
  pivot_wider(names_from = ONT, values_from = Freq, names_prefix = "Nanopore_") |>
  select(Sample, EM_Bin = EM, Nanopore_Low, Nanopore_Medium, Nanopore_High) |>
  arrange(Sample, EM_Bin)

# Print to console
cat("\n--- Formatted Confusion Matrices ---\n")
print(supp_table_df)

# 3. Save as a clean CSV for your Supplementary Materials Word doc
write_csv(supp_table_df, file.path(out_dir, "Supplementary_Table_S2_EM_Confusion_Matrices.csv"))

# Optional: Save raw R object in case you need it later
saveRDS(conf_mat_list, file.path(out_dir, "confusion_matrices_EM.rds"))

message("Processing complete! Summary and Confusion Matrices saved to: ", out_dir)