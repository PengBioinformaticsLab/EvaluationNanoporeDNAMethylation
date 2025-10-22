# ==============================================================
# 03_stats_Nano_EM.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script computes statistics comparing Nanopore and EM-Seq
# methylation datasets. It calculates:
# - Number of detected CpG sites at different coverage thresholds.
# - Overlap of CpG sites between both datasets.
# - Pearson correlation of methylation beta values.
#
# Dependencies:
# - Requires input data frames from Nanopore and EM-Seq (.rds files).
#
# Outputs:
# - Data frame containing CpG detection counts, overlaps, and correlations.
#
# ==============================================================


# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Preprocess Nanopore methylation data
#'
#' Computes beta values and generates a unique position identifier
#' for each CpG site.
#'
#' @param nano_df Data frame containing Nanopore CpG methylation data.
#' @return Data frame with added `beta` and `pos` columns.
preproc_nano <- function(nano_df){
  # Check for required columns
  required_cols <- c("NMod", "NOtherMod", "score", "chr", "start")
  missing_cols <- setdiff(required_cols, colnames(nano_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse=", ")))
  }
  
  # Convert relevant columns to numeric
  nano_df$NMod <- as.numeric(nano_df$NMod)
  nano_df$NOtherMod <- as.numeric(nano_df$NOtherMod)
  nano_df$score <- as.numeric(nano_df$score)
  
  # Handle division by zero
  nano_df$beta <- ifelse(nano_df$score == 0, NA, (nano_df$NMod + nano_df$NOtherMod) / nano_df$score)
  
  # Create unique position identifier
  nano_df$pos <- paste(nano_df$chr, (nano_df$start + 1), sep="_")
  
  return(nano_df)
}


#' Preprocess EM-Seq methylation data
#'
#' Computes beta values and generates a unique position identifier
#' for each CpG site.
#'
#' @param em_df Data frame containing EM-Seq CpG methylation data.
#' @return Data frame with added `beta` and `pos` columns.
preproc_em <- function(em_df) {
  # Check for required columns
  required_cols <- c("methy", "score", "chr", "start")
  missing_cols <- setdiff(required_cols, colnames(em_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse=", ")))
  }
  
  # Convert relevant columns to numeric
  em_df$methy <- as.numeric(em_df$methy)
  em_df$score <- as.numeric(em_df$score)
  
  # Handle division by zero
  em_df$beta <- ifelse(em_df$score == 0, NA, em_df$methy / em_df$score)
  
  # Create unique position identifier
  em_df$pos <- paste(em_df$chr, em_df$start, sep="_")
  
  return(em_df)
}


#' Compute CpG site statistics across coverage thresholds
#'
#' This function calculates:
#' - The number of detected CpG sites for both Nanopore and EM-Seq.
#' - Overlapping CpG sites between both datasets.
#' - Pearson correlation coefficient for shared CpG methylation values.
#'
#' @param nano_df Data frame with Nanopore CpG methylation data.
#' @param em_df Data frame with EM-Seq CpG methylation data.
#' @return A data frame containing CpG counts, overlaps, and correlation values.
count_cpgs <- function(nano_df, em_df) {
  
  # Check for required columns
  required_cols <- c("chr", "start", "score", "beta")
  if (!all(required_cols %in% colnames(nano_df))) {
    stop("nano_df is missing required columns: ", paste(setdiff(required_cols, colnames(nano_df)), collapse=", "))
  }
  if (!all(required_cols %in% colnames(em_df))) {
    stop("em_df is missing required columns: ", paste(setdiff(required_cols, colnames(em_df)), collapse=", "))
  }
  
  # Preprocess input datasets
  nano_df <- preproc_nano(nano_df)
  em_df <- preproc_em(em_df)
  
  # Initialize storage vectors
  nano_cpgs_detected <- NULL
  em_cpgs_detected <- NULL
  overlap_cpgs_detected <- NULL
  cor_df <- NULL
  
  # Define coverage thresholds
  cutoff_coverage <- c(1, 5, 10, 15, 20)
  
  for (cutoff in cutoff_coverage) {
    # Apply coverage threshold filtering
    nano_df <- nano_df[nano_df$score >= cutoff, ]
    em_df <- em_df[em_df$score >= cutoff, ]
    
    # Count CpG sites
    nano_cpgs_detected <- c(nano_cpgs_detected, nrow(nano_df))
    em_cpgs_detected <- c(em_cpgs_detected, nrow(em_df))
    
    # Identify overlapping CpG sites
    overlap_cpgs_detected <- c(overlap_cpgs_detected, sum(em_df$pos %in% nano_df$pos))
    
    # Extract shared CpG sites
    shared_pos <- em_df$pos[em_df$pos %in% nano_df$pos]
    
    data_em_sel <- em_df[em_df$pos %in% shared_pos, ]
    data_nano_sel <- nano_df[nano_df$pos %in% shared_pos, ]
    
    # Ensure order consistency
    if (sum(data_em_sel$pos == data_nano_sel$pos) != nrow(data_em_sel)) {
      stop("Error: Sample mismatch in CpG site selection.")
    }
    
    # Compute Pearson correlation for shared CpG sites
    cor_df <- c(cor_df, cor(data_em_sel$beta, data_nano_sel$beta, method = "pearson"))
  }
  
  # Compile results into a data frame
  cpgs_detected <- data.frame(
    coverage_threshold = paste0(cutoff_coverage, "X"),
    nano_cpgs_detected = nano_cpgs_detected,
    em_cpgs_detected = em_cpgs_detected,
    overlap_cpgs_detected = overlap_cpgs_detected,
    correlation = round(cor_df, 3)
  )
  
  return(cpgs_detected)
}

# ==============================================================
# UNIT TESTING
# ==============================================================


# Create synthetic test data
nano_test <- data.frame(
  chr = c("chr1", "chr1", "chr2"),
  start = c(100, 200, 300),
  NMod = c(5, 10, 15),
  NOtherMod = c(2, 3, 4),
  score = c(10, 20, 30)
)

em_test <- data.frame(
  chr = c("chr1", "chr1", "chr2"),
  start = c(101, 200, 300),
  methy = c(7, 13, 19),
  score = c(10, 20, 30)
)

# Test preprocessing functions
test_that("preproc_nano computes correct beta and pos", {
  result <- preproc_nano(nano_test)
  expect_equal(result$beta, c(0.7, 0.65, 0.6333), tolerance = 1e-3)
  expect_equal(result$pos, c("chr1_101", "chr1_201", "chr2_301"))
})

test_that("preproc_em computes correct beta and pos", {
  result <- preproc_em(em_test)
  expect_equal(result$beta, c(0.7, 0.65, 0.6333), tolerance = 1e-3)
  expect_equal(result$pos, c("chr1_101", "chr1_200", "chr2_300"))
})

# Test CpG counting function
test_that("count_cpgs correctly counts CpG sites and overlaps", {
  result <- count_cpgs(preproc_nano(nano_test),
                       preproc_em(em_test))
  
  expect_equal(result$nano_cpgs_detected[1], 3)
  expect_equal(result$em_cpgs_detected[1], 3)
  expect_equal(result$overlap_cpgs_detected[1], 1)  # Overlaps at chr1_200 and chr2_300
  
})





# ==============================================================
# GENERATE STATISTICS TABLE
# ==============================================================

#' Create a summary statistics table for CpG sites across multiple samples
#'
#' Computes CpG statistics for four sample pairs and writes to a CSV file.
#'
#' @param nano_1, em_1 Data frames for sample 1 (Nanopore, EM-Seq).
#' @param nano_2, em_2 Data frames for sample 2 (Nanopore, EM-Seq).
#' @param nano_3, em_3 Data frames for sample 3 (Nanopore, EM-Seq).
#' @param nano_4, em_4 Data frames for sample 4 (Nanopore, EM-Seq).
#' @return Saves a CSV file with CpG detection statistics.
create_stats_table <- function(nano_1, em_1,
                               nano_2, em_2,
                               nano_3, em_3,
                               nano_4, em_4) {
  
  # Compute CpG statistics for each sample
  stats_table <- rbind(
    count_cpgs(nano_1, em_1),
    count_cpgs(nano_2, em_2),
    count_cpgs(nano_3, em_3),
    count_cpgs(nano_4, em_4)
  )
  
  # Generate sample labels
  sample_name <- rep(paste0("Sample ", 1:4), each = 5)
  
  # Format the results
  stats_table <- cbind(sample_name, stats_table)
  stats_table <- data.frame(stats_table)
  
  # Define column headers
  colnames(stats_table) <- c(
    "Sample",
    "Coverage Threshold",
    "CpGs Detected (Nanopore)",
    "CpGs Detected (EM-Seq)",
    "CpGs Shared (Nanopore & EM-Seq)",
    "Correlation (Pearson)"
  )
  
  # Save results to CSV
  write.csv(stats_table, "output/tables/cpg_stats.csv", row.names = FALSE)
}


# ==============================================================
# EXECUTE STATISTICS COMPUTATION
# ==============================================================

# Compute statistics for all sample pairs
create_stats_table(nano_1, em_1,
                   nano_2, em_2,
                   nano_3, em_3,
                   nano_4, em_4)







