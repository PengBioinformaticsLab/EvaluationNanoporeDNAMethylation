#' Downsample Nanopore Methylation Data
#' 
#' This script takes a BED-methyl formatted RDS file from Nanopore sequencing
#' and downsamples it to a specified average coverage using a binomial distribution.
#' The function preserves the relative proportion of modified and canonical bases.
#' 
#' @param input_path Character string specifying the path to the .rds file containing the methylation data.
#' @param final_avg_cov Numeric value specifying the target average coverage after downsampling.
#' 
#' @return Saves the downsampled data as an .rds file in the same directory.
#' "'

# Downsampler function
# This function reads a BED-methyl formatted RDS file, computes the necessary downsampling ratio,
# and applies binomial sampling to reduce the coverage while maintaining relative proportions.
# 
# @param input_path Character, path to the input .rds file.
# @param final_cov Numeric, target mean coverage.
# @return Saves the downsampled dataframe to an .rds file.
downsample <- function(input_path, final_cov) {
  
  nano_df <- readRDS(file = input_path)
  
  # Find the downsampling ratio first: mean(nano_df$score) * X = final_cov |  X = ?
  ratio <- final_cov / mean(nano_df$score)
  
  # Downsample the dataframe by the ratio
  nano_df$NMod <- rbinom(
    rep(1, nrow(nano_df)),
    nano_df$NMod, ratio
  )
  nano_df$NCanonical <- rbinom(
    rep(1, nrow(nano_df)),
    nano_df$NCanonical, ratio
  )
  nano_df$NOtherMod <- rbinom(
    rep(1, nrow(nano_df)),
    nano_df$NOtherMod, ratio
  )
  nano_df$score <- nano_df$NMod + nano_df$NCanonical + nano_df$NOtherMod
  
  nano_df <- nano_df[nano_df$score > 0, ]
  
  return(nano_df)
}