# ==============================================================
# 04_modtype_stats.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script defines an S4 class (`NanoObject`) for handling 
# Nanopore CpG methylation data. It stores sample metadata, 
# CpG sites, and their overlap with genomic annotations such as
# CpG islands and promoters. It then summarizes the methylation data 
# based on genomic region.
#
# Dependencies:
# - Requires `GenomicRanges` for genomic region processing.
# - Assumes `cpg_annotation` is loaded in the environment.
#
# Outputs:
# - Summary table 
# - A summary table containing
#
# ==============================================================


# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Create an S4 object to store CpG methylation data
#'
#' This function defines and populates an S4 object (`NanoObject`) 
#' that stores CpG methylation information, sample metadata, and 
#' CpG site overlaps with genomic annotations.
#'
#' @param nano_df Data frame containing CpG methylation data.
#' @param sample_label Character string representing the sample ID.
#' @param tissue Character string indicating tissue type.
#' @return An S4 object (`NanoObject`) containing CpG site data and annotations.
make_s4 <- function(nano_df, sample_label, tissue) {
  
  nano_df$beta <- (nano_df$NMod + nano_df$NOtherMod) / nano_df$score
  
  setClass(
    "NanoObject",
    slots = list(
      sample = "character",
      tissue = "character",
      cpgs = "data.frame",
      island_cpgs_indices = "integer",
      promoter_cpgs_indices = "integer"
    )
  )
  
  mys4 <- new("NanoObject")
  
  
  mys4@sample <- sample_label
  mys4@tissue <- tissue
  
  mys4@cpgs <- nano_df
  
  
  region_indices <- function(nano_df, type) {
    nano_gr <- GRanges(
      seqnames = nano_df$chr,
      ranges = IRanges(nano_df$start, nano_df$end),
      NMod = nano_df$NMod,
      NOtherMod = nano_df$NOtherMod,
      score = nano_df$score,
      beta = nano_df$beta
    )
    
    # island_overlaps <- findOverlaps(nano_gr, cpg_annotation[cpg_annotation$type == "CpG_Island"])
    overlaps <- findOverlaps(nano_gr, cpg_annotation[cpg_annotation$type == type])
    
    return(queryHits(overlaps))
  }
  
  # "chr"        "start"      "end"        "modeCode"   "score"      "NMod"       "NCanonical" "NOtherMod"
  
  mys4@island_cpgs_indices <- region_indices(mys4@cpgs, "CpG_Island")
  mys4@promoter_cpgs_indices <- region_indices(mys4@cpgs, "Promoter")
  
  
  # Suggested usage:
  # nano_1@cpgs[queryHits(mys4@promoter_cpgs_indices ),]
  
  return(mys4)
}

# ================
#   Unit Testing
# ================

# Mock CpG methylation data
df <- data.frame(
  chr = c("chr1", "chr1", "chr2"),
  start = c(100, 200, 300),
  end = c(150, 250, 350),
  NMod = c(5, 10, 15),
  NOtherMod = c(2, 4, 6),
  score = c(10, 20, 30)
)

sample_label <- "Sample_X"
tissue <- "Liver"

# Create NanoObject instance
nano_obj <- make_s4(df, sample_label, tissue)

# Test initialization
test_that("NanoObject is initialized correctly", {
  expect_s4_class(nano_obj, "NanoObject")
  expect_equal(nano_obj@sample, sample_label)
  expect_equal(nano_obj@tissue, tissue)
  expect_true(is.data.frame(nano_obj@cpgs))
  expect_equal(nrow(nano_obj@cpgs), nrow(df))
  expect_type(nano_obj@island_cpgs_indices, "integer")
  expect_type(nano_obj@promoter_cpgs_indices, "integer")
})



#' Summarize methylation statistics for a given Nanopore dataset
#'
#' This function calculates key summary statistics for methylation data, including:
#' - The number of CpG sites.
#' - Mean, median, and standard deviation of methylation beta values.
#' - Mean coverage of CpG sites.
#'
#' @param nano_s4 An S4 object of class `NanoObject` containing CpG methylation data
#' @return A data frame with summarized methylation statistics.
methyl_summary <- function(nano_s4) {
  
  
  stats_helper <- function(cpgs_df){
    return(data.frame(
      "Number of 5mC Detected" = sum(cpgs_df[,"NMod"]),
      "Mean Beta 5mC" =  round(mean(cpgs_df[, "NMod"] / cpgs_df[, "score"]), 2),
      "Std Dev Beta 5mC" =  round(sd(cpgs_df[, "NMod"] / cpgs_df[, "score"]), 2),
      "Number of 5hmC Detected" = sum(cpgs_df[,"NOtherMod"]),
      "Mean Beta 5hmC" =  round(mean(cpgs_df[, "NOtherMod"] / cpgs_df[, "score"]), 2),
      "Std Dev Beta 5hmC" =  round(sd(cpgs_df[, "NOtherMod"] / cpgs_df[, "score"]), 2)
    ))
  }
  
  
  
  mod_summary <- rbind(
    stats_helper(nano_s4@cpgs),
    stats_helper(nano_s4@cpgs[nano_s4@island_cpgs_indices, ]),
    stats_helper(nano_s4@cpgs[nano_s4@promoter_cpgs_indices, ])
  )
  
  
  Sample <- rep(nano_s4@sample, 3)
  Tissue <-  rep(nano_s4@tissue, 3)
  Region <-  c("Global", "CpG Islands", "Promoter")
  
  
  mod_summary <- cbind(Sample, Tissue, Region, mod_summary)
  
  
  
  return(mod_summary)
}




# Test that the summary table has the correct format
test_that("methyl_summary returns a correctly formatted table", {
  # Create a mock NanoObject (assuming a constructor function exists)
  
  # Define the expected column names for the summary table
  expected_columns <- c("Sample", "Tissue", "Region", 
                        "Number.of.5mC.Detected", "Mean.Beta.5mC",
                        "Std.Dev.Beta.5mC", "Number.of.5hmC.Detected",
                        "Mean.Beta.5hmC")
  

  
  # Run the summary function
  summary_table <- methyl_summary(nano_obj)
  
  # Check if the output is a data frame
  expect_s3_class(summary_table, "data.frame")
  
  # Check if all expected columns are present
  expect_true(all(expected_columns %in% colnames(summary_table)))
  
  # Ensure the table is not empty
  expect_gt(nrow(summary_table), 0)
})



# ==============================================================
# EXECUTE S4 OBJECT CREATION / SUMMARIZATION
# ==============================================================

sample_report <- function(df, sample_label, tissue) {
  nano_s <- make_s4(df, sample_label, tissue)
  return(methyl_summary(nano_s))
}

mod_summary <- rbind(sample_report(nano_1, "1", "Brain"),
                     sample_report(nano_2, "2", "Brain"),
                     sample_report(nano_3, "3", "Brain"),
                     sample_report(nano_4, "4", "Brain"),
                     sample_report(nano_a, "A", "Blood"),
                     sample_report(nano_b, "B", "Blood"),
                     sample_report(nano_c, "C", "Blood"),
                     sample_report(nano_d, "D", "Blood"))


write.csv(mod_summary, file = "output/tables/5mC_5hmC_table.csv")

