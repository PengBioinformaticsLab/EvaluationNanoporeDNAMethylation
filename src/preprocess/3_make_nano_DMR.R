# ==============================================================
# 3_make_nano_DMR.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script processes haplotype-specific Nanopore `.bedmethyl.gz` files,
# applies coverage filtering, and calls differentially methylated regions (DMRs)
# using DSS. The results are saved as `.rds` files for downstream analysis.
#
# Dependencies:
# - Requires `DSS` and `bsseq` for DMR calling.
# - Requires `.bedmethyl.gz` files generated from Epi2Me pipeline.
#
# Outputs:
# - Processed haplotype `.rds` files (hap1.rds, hap2.rds).
# - DMR results saved as `dmrs.rds`.
#
# ==============================================================



# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Load and preprocess a haplotype-specific `.bedmethyl.gz` file
#'
#' This function reads a `.bedmethyl.gz` file, filters for CpG sites,
#' applies a minimum coverage threshold, and converts it into a BSseq object.
#'
#' @param bedmethyl_file_hap Character. Path to the haplotype `.bedmethyl.gz` file.
#' @param min_cov Numeric. Minimum coverage threshold for filtering CpG sites.
#' @return A BSseq object containing filtered methylation data.
loadHap <- function(bedmethyl_file_hap, min_cov = 1) {
  cn <- c(
    "chr", "start", "end", "modiCode", "score", "strand", "start2", "end2",
    "color", "NVCov", "percentModi", "NMod", "NCanonical", "NOtherMod",
    "NDelete", "NFail", "NDiff", "NNoCall"
  )
  
  cn_sel <- c(
    "chr", "start", "end", "modiCode", "score", "strand", "percentModi",
    "NMod", "NCanonical", "NOtherMod"
  )
  
  bedmethy_hap <- read.table(
    
    gzfile(bedmethyl_file_hap),
    header = FALSE, sep = "", col.names = cn
  )
  
  bedmethy_hap <- bedmethy_hap[, cn_sel]
  
  # standard chromosomes
  bedmethy_hap <- bedmethy_hap[bedmethy_hap$chr %in% paste0("chr", c(1:22, "X", "Y")), ]
  # 5mC as reference modcode
  bedmethy_hap <- bedmethy_hap[bedmethy_hap$modiCode == "m", ]
  
  # threshold on coverage
  bedmethy_hap <- bedmethy_hap[bedmethy_hap$score >= min_cov, ]
  
  # change to a BSobj from BSS package to use the DSS package
  bedmethy_hap <- transform(bedmethy_hap, chr = chr, pos = start, N = score, X = NMod + NOtherMod)[, c("chr", "pos", "N", "X")]
  
  return(bedmethy_hap)
}

#' Call differentially methylated regions (DMRs) between haplotypes
#'
#' This function performs DMR analysis using DSS on two haplotype-specific 
#' BSseq objects, applying statistical filtering.
#'
#' @param BS_hap1 BSseq object for haplotype 1.
#' @param BS_hap2 BSseq object for haplotype 2.
#' @return Data frame containing identified DMRs.
call_dmrs <- function(hap1, hap2) {
  BS_obj <- makeBSseqData(
    list(
      hap1,
      hap2
    ),
    c("hap1", "hap2")
  )
  
  ### The below warning is okay and can be ignored after you sort
  # Warning: CG positions in chromosome chr1 is not ordered. Reorder CG sites.
  
  BS_obj <- orderBSseq(BS_obj, seqOrder = c(paste0("chr", 1:22), "chrX", "chrY"))
  
  # The ncores = 4 call below is important. The supercomputer will kick us off otherwise.
  # also this might take forever.
  dmlTest.sm <- DMLtest(BS_obj, group1 = c("hap1"), group2 = c("hap2"), equal.disp = FALSE, smoothing = TRUE, smoothing.span = 500, ncores = 4)
  
  dmrs <- callDMR(dmlTest.sm) # Default p-value threshold is 1e-5
  
  return(dmrs)
}


# ==============================================================
# MAIN FUNCTION TO PROCESS HAPLOTYPES AND CALL DMRs
# ==============================================================

#' Process haplotype `.bedmethyl.gz` files and call DMRs
#'
#' This function reads haplotype-specific methylation data, applies filtering,
#' calls DMRs, and saves the results to `.rds` files.
#'
#' @param sample_input_dir Character. Directory containing `.bedmethyl.gz` files.
#' @param out_dir Character. Directory for output `.rds` files.
#' @param samplename Character. Base name for output files.
#' @param min_cov Numeric. Minimum coverage threshold for filtering.
#' @return Saves processed haplotype `.rds` files and `dmrs.rds` to disk.
dmrs_pipeline <- function(sample_input_dir, out_dir, samplename, min_cov = 5){
  
  ### Use the input directory to get the hap1 and hap2 files...
  hap1_file <- file.path(sample_input_dir, get_bedmethyl_files(sample_input_dir)['h1'])
  hap2_file <- file.path(sample_input_dir, get_bedmethyl_files(sample_input_dir)['h2'])
  
  
  ### hap1
  hap1 <- loadHap(hap1_file, min_cov = min_cov)
  ### hap2
  hap2 <- loadHap(hap2_file, min_cov = min_cov)
  
  ### calculate DMRS
  dmrs <- call_dmrs(hap1 = hap1, hap2 = hap2)
  
  
  saveRDS(object = dmrs, file = file.path(out_dir, paste0("dmrs_", min_cov, "X_", samplename, ".rds")))
  saveRDS(object = hap1, file = file.path(out_dir, paste0("hap1_", min_cov,"X_", samplename, ".rds")))
  saveRDS(object = hap2, file = file.path(out_dir, paste0("hap2_", min_cov,"X_", samplename, ".rds")))
}

# ==============================================================
# EXECUTE DMR ANALYSIS (EXAMPLE USAGE)
# ==============================================================

# Example usage:
# dmrs_pipeline(sample_input_dir = sample_51_dir, out_dir = dmrs_out_dir, samplename = "1", min_cov = 5)






