# ==============================================================
# 2_save_EM_rds.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script processes Bismark output files (.cov.gz) for EM-Seq data.
# It filters for CpG sites, merges strands, and saves the processed 
# data as an .rds file.
#
# Dependencies:
# - Requires `data.table` for fast data reading.
# - Requires `RDS` format for efficient data storage.
#
# Outputs:
# - Processed EM-Seq data saved as an .rds file.
#
# ==============================================================


# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Read and preprocess an EM-Seq file
#'
#' This function reads an EM-Seq methylation file (.cov.gz) using `data.table`,
#' filters for CpG sites, and removes unnecessary columns.
#'
#' @param em_filepath Character. File path to the EM-Seq sample.
#' @return Data table with processed EM-Seq methylation data.
read_EM <- function(em_filepath){
  
  #column names
  cn <- c("chr",
          "start",
          "strand",
          "methy",
          "unmethy",
          "C_context") #CpG or alternative context)
  
  #EM datatable object. 7th column doesnt matter
  em_dt <- data.table::fread(em_filepath, drop = 7, col.names = cn)
  
  #Filter for CpGs using datatable methods, then drop the column to save precious memory
  em_dt <- em_dt[C_context == "CG", ]
  em_dt <- em_dt[, C_context := NULL]
  
  return(em_dt)
}


# ==============================================================
# COMBINE STRANDS FUNCTION
# ==============================================================

#' Merge complementary DNA strands in EM-Seq data
#'
#' This function processes a data frame of EM-Seq methylation data by:
#' - Adjusting coordinates for the negative strand (`-1` to match CpG site location).
#' - Combining methylation and unmethylation counts for both strands.
#'
#' @param em_df Data frame containing processed EM-Seq methylation data.
#' @return Data frame with merged strand information.
strand_combo <- function(em_df) {
  
  em_final <- NULL
  
  std_chr <- paste0("chr", c(1:22,'X','Y'))
  
  for (c in std_chr) {
    # Combine a chromosome at a time
    print(c)
    em_sel <- em_df[em_df$chr == c, c("start", "strand", "methy", "unmethy")]
    
    # Subtract one from the start position for the '-' strand.
    em_sel[em_sel$strand == "-", ]$start <- em_sel[em_sel$strand == "-", ]$start - 1
    
    # Drop strand column
    em_sel$strand <- NULL
    
    # Sum methy, unmethy based on combined strand shared start positions
    combined_methy <- tapply(em_sel$methy, em_sel$start, sum)
    combined_unmethy <- tapply(em_sel$unmethy, em_sel$start, sum)
    
    #### Coerce the above methylation counts into a df to merge each on 'start' position
    # methy df
    combined_methy <- as.data.frame(as.table(combined_methy), stringsAsFactors = FALSE)
    names(combined_methy) <- c("start", "methy")
    
    # unmethy df
    combined_unmethy <- as.data.frame(as.table(combined_unmethy), stringsAsFactors = FALSE)
    names(combined_unmethy) <- c("start", "unmethy")
    
    # Merge both
    em_sel <- merge(combined_methy, combined_unmethy, by = "start")

    
    em_sel <- transform(
      em_sel,
      chr = c,
      start = as.integer(start), 
      score = methy + unmethy
    )
    
    # Select only the desired columns
    em_sel <- em_sel[c("chr", "start", "score", "methy", "unmethy")]
    
    # Arrange the rows by the 'start' column
    em_sel <- em_sel[order(em_sel$start), ]
    
    # Combine with rest of processed chromosome sections
    em_final <- rbind(em_final, em_sel)
    
   # Remove the processed chr from the original EM_df.
    #lower memory usage
    em_df <- em_df[em_df$chr != c,]
    gc()
  }
  
  em_final$beta <- round((em_final$methy / em_final$score),3)
  
  
  return(em_final)
} 

# SAVE PROCESSED DATA
# ==============================================================

#' Save processed EM-Seq data as an RDS file
#'
#' This function reads an EM-Seq `.cov.gz` file, processes it, 
#' merges strands, and saves the final dataset as an RDS file.
#'
#' @param em_filepath Character. Path to the EM-Seq file.
#' @param output_filepath Character. Path for saving the .rds file.
#' @return Saves processed data to `output_filepath`.
em_proc <- function(em_filepath, output_filepath){
  
  em <- read_EM(em_filepath)
  em <- strand_combo(em)
  
  saveRDS(em, file = output_filepath)
}



# ==============================================================
# EXECUTE PIPELINE (EXAMPLE USAGE)
# ==============================================================

# Example usage:
# save_em_rds("data/EM_sample.cov.gz", "output/em_sample.rds")






