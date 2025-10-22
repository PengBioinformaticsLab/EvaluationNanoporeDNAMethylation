# ==============================================================
# 01_filepaths.R
# Author: Steven Brooks
# Date: 02/15/2025
# --------------------------------------------------------------
# Description:
# This script loads methylation data stored in `.rds` files for 
# both Nanopore and EM-Seq datasets. It also constructs genomic 
# ranges (`GRanges`) for downstream analysis. The CpG annotation 
# file is loaded, and overlaps with CpG islands are determined.
#
# Dependencies:
# - Requires GenomicRanges and IRanges for genomic operations.
# - `.rds` files must be available in the specified directories.
#
# Outputs:
# - GRanges object (`nano_gr`) for methylation analysis.
# - CpG annotation (`cpg_annotation`) for genomic overlap.
#
# ==============================================================

# Function to check and load .rds files safely
load_rds_safe <- function(filepath) {
  if (!file.exists(filepath)) {
    warning(paste("Warning: Missing .rds file -", filepath))
    return(NULL)
  }
  
  # Attempt to read the .rds file and catch corruption errors
  result <- tryCatch(
    {
      readRDS(filepath)
    },
    error = function(e) {
      stop(paste("Error: Corrupt .rds file or invalid format -", filepath, "-", e$message))
    }
  )
  
  return(result)
}



# ==============================================================
# LOADING .rds FILES
# ==============================================================


# Nanopore methylation datasets
nano_1 <- load_rds_safe(file = 'data/nanopore/nano_51.rds')  # Sample 1
nano_2 <- load_rds_safe(file = 'data/nanopore/nano_69.rds')  # Sample 2
nano_3 <- load_rds_safe(file = 'data/nanopore/nano_168.rds') # Sample 3
nano_4 <- load_rds_safe(file = 'data/nanopore/nano_266.rds') # Sample 4

# EM-Seq methylation datasets
em_1 <- load_rds_safe(file = 'data/em-seq/em_51.rds')   # Sample 1
em_2 <- load_rds_safe(file = 'data/em-seq/em_69.rds')   # Sample 2
em_3 <- load_rds_safe(file = 'data/em-seq/em_168.rds')  # Sample 3
em_4 <- load_rds_safe(file = 'data/em-seq/em_266.rds')  # Sample 4

# Additional Nanopore data from EPIC comparison
nano_a <- load_rds_safe('../nanopore_epic/data/sample_4059246541_analysis/s4059246541_combined.rds')
nano_b <- load_rds_safe('../nanopore_epic/data/sample_4059249500_analysis/s4059249500_combined.rds')
nano_c <- load_rds_safe('../nanopore_epic/data/sample_4059246577_analysis/s4059246577_combined.rds')
nano_d <- load_rds_safe('../nanopore_epic/data/sample_4059249475_Analysis/s4059249475_combined.rds')

# ==============================================================
# CONVERTING DATA TO GENOMIC RANGES (GRanges)
# ==============================================================

# Construct a GRanges object for Nanopore methylation data
# This allows for efficient genomic range operations.
nano_gr <- GRanges(seqnames = nano_a$chr,
                   ranges = IRanges(nano_a$start, nano_a$end),
                   NMod = nano_a$NMod,             # Number of modified CpGs
                   NOtherMod = nano_a$NOtherMod,   # Number of other modifications
                   beta = (nano_a$NMod + nano_a$NOtherMod) / nano_a$score, # Methylation ratio
                   score = nano_a$score)           # Coverage score

# ==============================================================
# LOADING REFERENCE DATA
# ==============================================================

# Load CpG annotation file for genomic feature mapping
cpg_annotation <- load_rds_safe(file = 'data/reference/cpg_annotation.rds')

# ==============================================================
# IDENTIFYING OVERLAPS WITH CpG ISLANDS
# ==============================================================

# Find overlaps between Nanopore methylation data and CpG Islands
island_overlaps <- findOverlaps(nano_gr, cpg_annotation[cpg_annotation$type == "CpG_Island"])

# ==============================================================
# TESTING
# ==============================================================


# Function to validate GRanges object
validate_granges <- function(gr_obj, name) {
  if (!inherits(gr_obj, "GRanges")) {
    stop(paste("Error:", name, "is not a valid GRanges object."))
  }
  
  message(paste("Validation passed for", name))
}

# Function to validate CpG annotation overlap
validate_cpg_overlap <- function(meth_data, cpg_annotation) {
  if (!inherits(meth_data, "GRanges")) {
    stop("Error: Methylation dataset is not a valid GRanges object.")
  }
  if (!inherits(cpg_annotation, "GRanges")) {
    stop("Error: CpG annotation is not a valid GRanges object.")
  }
  
  overlaps <- findOverlaps(meth_data, cpg_annotation)
  if (length(overlaps) == 0) {
    warning("Warning: No overlaps found between methylation data and CpG annotation.")
  } else {
    message("CpG annotation overlap validation passed. Overlaps detected.")
  }
}

# Ensure critical files are loaded
if (is.null(nano_1) || is.null(nano_2) || is.null(cpg_annotation)) {
  stop("Error: One or more critical .rds files are missing. Please check file paths.")
}

# Continue with processing if all critical files are available
message("The .rds files are successfully loaded.")


# Validate GRanges objects
validate_granges(nano_gr, "nano_a")
validate_granges(cpg_annotation, "cpg_annotation")

# Validate CpG annotation overlap
validate_cpg_overlap(nano_gr, cpg_annotation)
