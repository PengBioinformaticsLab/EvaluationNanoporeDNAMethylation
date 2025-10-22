# ==============================================================
# 10_DMR_stats.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script reads differentially methylated region (DMR) data from multiple 
# samples, summarizes key statistics (count, mean length, mean nCG count),
# and exports a summary table.
#
# Dependencies:
# - Requires `dplyr` for data manipulation.
# - Requires `.rds` files containing DMR datasets.
#
# Outputs:
# - A CSV summary file with statistics for each sample.
#
# ==============================================================


# ==============================================================
# DEFINE FILE PATHS AND SAMPLE INFORMATION
# ==============================================================

# Sample identifiers
samples <- c("1", "2", "3", "4", "A", "B", "C", "D")

# Sample sources (e.g., Brain vs. Blood)
sources <- c(rep("Brain", 4), rep("Blood", 4))

# Generate file paths to DMR datasets
file_paths <- paste0("data/dmrs/dmrs_5X_", samples, ".rds")

# ==============================================================
# LOAD DMR DATASETS
# ==============================================================

#' Read and load DMR datasets
#'
#' This function reads a list of `.rds` files containing DMR information.
#'
#' @param file_paths Character vector containing paths to `.rds` files.
#' @return A list of data frames, each containing DMR data for a sample.
read_dmrs <- function(file_paths) {
  lapply(file_paths, readRDS)
}

# Read all DMR datasets
dmrs_list <- read_dmrs(file_paths)


# ==============================================================
# SUMMARIZE DMR STATISTICS
# ==============================================================

#' Summarize key statistics for a given DMR dataset
#'
#' This function calculates the total count of DMRs, mean length, 
#' and mean number of CpG sites (nCG).
#'
#' @param dmrs Data frame containing DMRs.
#' @return Data frame with summary statistics (count, mean length, mean nCG).
summarize_dmrs <- function(dmrs) {
  dmrs %>%
    summarise(
      Count = n(),  # Total number of DMRs
      Mean_length = round(mean(length), 2),  # Mean DMR length
      Mean_nCG = round(mean(nCG), 2)  # Mean number of CpGs in DMRs
    )
}

# Apply summary function to all datasets and combine results
smry_table <- do.call(rbind, lapply(dmrs_list, summarize_dmrs))

# Add sample and source information to the summary table
smry_table <- cbind(sample = samples, source = sources, smry_table)


# ==============================================================
# EXPORT SUMMARY TABLE
# ==============================================================

# Print the summary table for review
print(smry_table)

# Save the summary table to a CSV file
write.csv(smry_table, file = "output/tables/DMRs_summary.csv", row.names = FALSE)