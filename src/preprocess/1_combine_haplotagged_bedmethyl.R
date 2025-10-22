# Note: Try modkit bedmethyl merge in the official docs.
#the below script was written before modkit released bedmethyl.gz merging

# ==============================================================
# 1_combine_haplotagged_bedmethyl.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script merges haplotype 1, haplotype 2, and ungrouped bedmethyl 
# data into a single file that contains combined 5mC and 5hmC methylation 
# information.

#
# Dependencies:
# - Requires bedmethyl files in `.gz` format.
# - Requires `dplyr` for data manipulation.
#
# Outputs:
# - A merged bedmethyl dataset.
#
# ==============================================================

#Function description:
#Merges haplotype 1, haplotype 2, and ungrouped bedmethyl data
#into a single file with combined 5mC and 5hmC information.

#' combine hap1, hap2, and ungrouped data into single bedmethyl file with 5mc and 5hmc combined

#'

#' @param bedmethyl_file_hap1 bedmethyl file from haplotype 1

#' @param bedmethyl_file_hap2 bedmethyl file from haplotype 2

#' @param bedmethyl_file_ungrouped ungrouped bedmethyl file

#' @param filename_out name of the file to output results. If it is NULL, return the resutls.
#' 


mergeHapsAllMethylation <- function(
    
  bedmethyl_file_hap1, bedmethyl_file_hap2,
  
  bedmethyl_file_ungrouped, filename_out = NULL
  
){

#Define column names for the input files, including chromosome, position, 
#modification details, and coverage information.
  
  cn <- c(
    
    "chr", "start", "end", "modiCode", "score", "strand", "start2", "end2",
    
    "color", "NVCov", "percentModi", "NMod", "NCanonical", "NOtherMod",
    
    "NDelete", "NFail", "NDiff", "NNoCall"
    
  )
  
  
#Select a subset of columns to retain in the final dataset.
  cn_sel <- c(
    
    "chr", "start", "end", "modiCode", "score", "strand", "percentModi",
    
    "NMod", "NCanonical", "NOtherMod"
    
  )
  
#Read the bedmethyl files for haplotype 1, haplotype 2, and ungrouped data from compressed .gz format.
  
  bedmethy_hap1 <- read.table(
    
    gzfile(bedmethyl_file_hap1),
    
    header = FALSE, sep = "", col.names = cn
    
  )  
  
  
  
  bedmethy_hap2 <- read.table(
    
    gzfile(bedmethyl_file_hap2),
    
    header = FALSE, sep = "", col.names = cn
    
  )
  
  
  
  bedmethy_ungrouped <- read.table(
    
    gzfile(bedmethyl_file_ungrouped),
    
    header = FALSE, sep = "", col.names = cn
    
  )
  
#Keep only the selected columns in each dataset.
  
  bedmethy_hap1 <- bedmethy_hap1[,cn_sel]
  
  bedmethy_hap2 <- bedmethy_hap2[,cn_sel]
  
  bedmethy_ungrouped <- bedmethy_ungrouped[,cn_sel]
  
  
#Debugging lines (commented out): Print modification type counts for each dataset.
  
  # print(table(bedmethy_hap1$modiCode))
  
  # print(table(bedmethy_hap2$modiCode))
  
  # print(table(bedmethy_ungrouped$modiCode))
  
#Define the chromosome list (excluding mitochondrial "MT" since its CpGs are only in the ungrouped file).
# not including "MT" here, CpGs on MT are all in ungrouped file
  
  chr_all <- paste0("chr", c(1:22, "X", "Y"))
  
#Initialize rlt_final as NULL to store final results.
  
  rlt_final <- NULL
  
  
#Loop over each chromosome and print the current chromosome being processed.
  for(chr in chr_all){
    
    print(chr)
    
    #Extract data for CpG sites (modification code "m") for haplotype 1, haplotype 2, and ungrouped data.
    
    h1 <- bedmethy_hap1[bedmethy_hap1$chr == chr & bedmethy_hap1$modiCode == "m", ]
    
    h2 <- bedmethy_hap2[bedmethy_hap2$chr == chr & bedmethy_hap2$modiCode == "m", ]
    
    ug <- bedmethy_ungrouped[bedmethy_ungrouped$chr == chr & bedmethy_ungrouped$modiCode == "m", ]
    
    
    
    start_all <- sort(unique(c(h1$start, h2$start, ug$start)))
    
    
    
    rlt <- data.frame(
      
      chr = chr,
      
      start = start_all,
      
      end = 0,
      
      modeCode = "m",
      
      score = 0,
      
      NMod = 0,
      
      NCanonical = 0,
      
      NOtherMod = 0,
      
      stringsAsFactors = FALSE
      
    )
    
    
    
    # add h1
    
    idx <- match(h1$start, rlt$start)
    
    rlt$end[idx] <- h1$end
    
    #rlt$score[idx] <- rlt$score[idx] + h1$score
    
    rlt$NMod[idx] <- rlt$NMod[idx] + h1$NMod
    
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + h1$NCanonical
    
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + h1$NOtherMod
    
    
    
    # add h2
    
    idx <- match(h2$start, rlt$start)
    
    rlt$end[idx] <- h2$end
    
    #rlt$score[idx] <- rlt$score[idx] + h2$score
    
    rlt$NMod[idx] <- rlt$NMod[idx] + h2$NMod
    
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + h2$NCanonical
    
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + h2$NOtherMod
    
    
    
    # add ug
    
    idx <- match(ug$start, rlt$start)
    
    rlt$end[idx] <- ug$end
    
    #rlt$score[idx] <- rlt$score[idx] + ug$score
    
    rlt$NMod[idx] <- rlt$NMod[idx] + ug$NMod
    
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + ug$NCanonical
    
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + ug$NOtherMod
    
    
    
    if(chr == "chr1"){
      
      rlt_final <- rlt
      
    } else {
      
      rlt_final <- rbind(rlt_final, rlt)
      
    }
    
  }
  
  
  #Side step any issues with unscored NFail NDiff, etc...
  rlt_final$score <- rlt_final$NMod + rlt_final$NCanonical + rlt_final$NOtherMod
  
  rlt_final$beta <- round((rlt_final$NMod + rlt_final$NOtherMod)/(rlt_final$score),3)
  
  if(is.null(filename_out)){
    
    return(rlt_final)
    
  }
  
  
  
  
  saveRDS(rlt_final, file = filename_out)
  
}



#'   #Given a sample directory, return a vector of file paths of hap1, hap2, ungrouped, and a sample_X.rds filepath
#'   Helper function for the above mergeHap
#' @param sample_dir directory containing hap1, hap2, and ungrouped .bedmethyl files
get_bedmethyl_files <- function(sample_dir){
  
  hap1_file <- list.files(path = sample_dir, pattern = "1\\.wf_mods\\.bedmethyl\\.gz|1\\.bedmethyl\\.gz")
  hap2_file <- list.files(path = sample_dir, pattern = "2\\.wf_mods\\.bedmethyl\\.gz|2\\.bedmethyl\\.gz")
  ungrouped_file  <- list.files(path = sample_dir, pattern = "ungrouped\\.wf_mods\\.bedmethyl\\.gz|ungrouped\\.bedmethyl\\.gz") 
  sample_name <- paste0('nano_', gsub(".*_(\\d+)_.*", "\\1", hap1_file), '.rds')
  
  
  return(c(h1 = hap1_file,h2 =  hap2_file,ungr =  ungrouped_file, sample_name =  sample_name))
}









