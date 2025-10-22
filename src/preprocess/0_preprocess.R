# ==============================================================
# 0_preprocess.R
# Author: Steven Brooks
# Date: 02/19/24
# --------------------------------------------------------------
# Description:
# This script sets up file paths and directories for processing 
# Nanopore and EM-Seq methylation data. It defines paths to raw data 
# and output directories to be used in downstream analysis.
#
# Dependencies:
# - Requires structured directories for Nanopore and EM-Seq files.
#
# Outputs:
# - Standardized paths for input and processed methylation data.
#
# ==============================================================

library(DSS)
library(bsseq)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(rtracklayer)

# For reproducibility
set.seed(123)


#### File paths needed:

# Define the base directory for Nanopore files
nano_files_path <- 'data/Dec2024_comparison/Nanopore'


# Define directories for specific Nanopore sample analyses
sample_51_dir <- file.path(nano_files_path, 'Su349_BTRC_51_Control_M_analysis')
sample_69_dir <- file.path(nano_files_path, 'Su395_BTRC_69_Control_M_analysis')
sample_168_dir <- file.path(nano_files_path, 'Su636_BTRC_168_Control_M_merge_analysis')
sample_266_dir <- file.path(nano_files_path, 'Su779_BTRC_266_AUD_M_analysis')



# Define the base directory for Epic Nanopore files
nano_files_path_epic <- '../nanopore_epic/data'

# Define directories for specific Epic Nanopore sample analyses
sample_541_dir <- file.path(nano_files_path_epic, 'sample_4059246541_analysis')
sample_577_dir <- file.path(nano_files_path_epic, 'sample_4059246577_analysis')
sample_475_dir <- file.path(nano_files_path_epic, 'sample_4059249475_Analysis')
sample_500_dir <- file.path(nano_files_path_epic, 'sample_4059249500_analysis')


### Output directory for processed nanopore data
nano_data_dir <- 'data/nanopore'

########## EM Pathing

# Define input directory for EM-Seq data
em_input_dir <- 'data/Dec2024_comparison/EM'

# Define paths for individual EM-Seq sample files
em_51_path <- file.path(em_input_dir,
                        'BTRC_51_BTRC_51_Control.CX_report.txt.gz')

em_69_path <- file.path(em_input_dir,
                        'BTRC_69_BTRC_69_Control.CX_report.txt.gz')

em_168_path <- file.path(em_input_dir,
                         'BTRC_168_BTRC_168_Control.CX_report.txt.gz')

em_266_path <- file.path(em_input_dir,
                         'BTRC_266_BTRC_266_AUD.CX_report.txt.gz')

### Output directory for EM-Seq processed data
em_data_dir <- 'data/em-seq'




# Directory for DMRs (Differentially Methylated Regions) analysis
dmrs_sample_dir <- "../nanopore_epic/data/sample_4059249500_analysis"
dmrs_out_dir <- "data/dmrs"

# Directory for reference genome or related files
ref_dir ='data/reference'



# ==============================================================
# FUNCTION TO VALIDATE FILE PATHS
# ==============================================================

#' Check if required data files exist
#'
#' This function checks whether the defined file paths for both 
#' Nanopore and EM-Seq datasets exist before proceeding.
#'
#' @param paths A character vector containing file paths to check.
#' @return Prints a message for each missing file.
validate_paths <- function(paths) {
  for (path in paths) {
    if (!file.exists(path)) {
      warning(paste("Warning: File not found:", path))
    }
  }
}

# Validate paths for Nanopore data
validate_paths(c(sample_51_dir, sample_69_dir, sample_168_dir, sample_266_dir,
                 sample_541_dir, sample_577_dir, sample_475_dir, sample_500_dir))

# Validate paths for EM-Seq data
validate_paths(c(em_51_path, em_69_path, em_168_path, em_266_path))



##### Load preprocessing script for merging haplotype-specific bedmethyl files
source('src/preprocess/1_combine_haplotagged_bedmethyl.R')

######### Process EM-Seq sample pairs by merging haplotype methylation data


##### Make .rds files:
source('src/preprocess/1_combine_haplotagged_bedmethyl.R')

######### EM-Seq sample pairs
# 51
mergeHapsAllMethylation(bedmethyl_file_hap1 = file.path(sample_51_dir, get_bedmethyl_files(sample_51_dir)['h1']),
                        bedmethyl_file_hap2 = file.path(sample_51_dir, get_bedmethyl_files(sample_51_dir)['h2']),
                        bedmethyl_file_ungrouped = file.path(sample_51_dir, get_bedmethyl_files(sample_51_dir)['ungr']),
                        filename_out = file.path(nano_data_dir,
                                                 get_bedmethyl_files(sample_51_dir)['sample_name']))

gc()

# 69
mergeHapsAllMethylation(bedmethyl_file_hap1 = file.path(sample_69_dir, get_bedmethyl_files(sample_69_dir)['h1']),
                        bedmethyl_file_hap2 = file.path(sample_69_dir, get_bedmethyl_files(sample_69_dir)['h2']),
                        bedmethyl_file_ungrouped = file.path(sample_69_dir, get_bedmethyl_files(sample_69_dir)['ungr']),
                        filename_out = file.path(nano_data_dir,
                                                 get_bedmethyl_files(sample_69_dir)['sample_name']))

gc()
# 168
mergeHapsAllMethylation(bedmethyl_file_hap1 = file.path(sample_168_dir, get_bedmethyl_files(sample_168_dir)['h1']),
                        bedmethyl_file_hap2 = file.path(sample_168_dir, get_bedmethyl_files(sample_168_dir)['h2']),
                        bedmethyl_file_ungrouped = file.path(sample_168_dir, get_bedmethyl_files(sample_168_dir)['ungr']),
                        filename_out = file.path(nano_data_dir,
                                                 get_bedmethyl_files(sample_168_dir)['sample_name']))

gc()

# 266
mergeHapsAllMethylation(bedmethyl_file_hap1 = file.path(sample_266_dir, get_bedmethyl_files(sample_266_dir)['h1']),
                        bedmethyl_file_hap2 = file.path(sample_266_dir, get_bedmethyl_files(sample_266_dir)['h2']),
                        bedmethyl_file_ungrouped = file.path(sample_266_dir, get_bedmethyl_files(sample_266_dir)['ungr']),
                        filename_out = file.path(nano_data_dir,
                                                 get_bedmethyl_files(sample_266_dir)['sample_name']))
gc()

######### EPIC sample pairs

# 475
mergeHapsAllMethylation(bedmethyl_file_hap1 = file.path(sample_475_dir, get_bedmethyl_files(sample_475_dir)['h1']),
                        bedmethyl_file_hap2 = file.path(sample_475_dir, get_bedmethyl_files(sample_475_dir)['h2']),
                        bedmethyl_file_ungrouped = file.path(sample_475_dir, get_bedmethyl_files(sample_475_dir)['ungr']),
                        filename_out = file.path(nano_data_dir, 'nano_475.rds'))

gc()

# 500
mergeHapsAllMethylation(bedmethyl_file_hap1 = file.path(sample_500_dir, get_bedmethyl_files(sample_500_dir)['h1']),
                        bedmethyl_file_hap2 = file.path(sample_500_dir, get_bedmethyl_files(sample_500_dir)['h2']),
                        bedmethyl_file_ungrouped = file.path(sample_500_dir, get_bedmethyl_files(sample_500_dir)['ungr']),
                        filename_out = file.path(nano_data_dir, 'nano_500.rds'))

gc()
# 541
mergeHapsAllMethylation(bedmethyl_file_hap1 = file.path(sample_541_dir, get_bedmethyl_files(sample_541_dir)['h1']),
                        bedmethyl_file_hap2 = file.path(sample_541_dir, get_bedmethyl_files(sample_541_dir)['h2']),
                        bedmethyl_file_ungrouped = file.path(sample_541_dir, get_bedmethyl_files(sample_541_dir)['ungr']),
                        filename_out = file.path(nano_data_dir, 'nano_541.rds'))

gc()
# 577
mergeHapsAllMethylation(bedmethyl_file_hap1 = file.path(sample_577_dir, get_bedmethyl_files(sample_577_dir)['h1']),
                        bedmethyl_file_hap2 = file.path(sample_577_dir, get_bedmethyl_files(sample_577_dir)['h2']),
                        bedmethyl_file_ungrouped = file.path(sample_577_dir, get_bedmethyl_files(sample_577_dir)['ungr']),
                        filename_out = file.path(nano_data_dir, 'nano_577.rds'))
gc()


# EM-Seq preprocess

source('src/preprocess/2_save_EM_rds.R')


# Process EM-Seq samples and save as .rds
em_proc(em_filepath = em_51_path, output_filepath = file.path(em_data_dir, 'em_51.rds'))
gc()
em_proc(em_filepath = em_69_path, output_filepath = file.path(em_data_dir, 'em_69.rds'))
gc()
em_proc(em_filepath = em_168_path, output_filepath = file.path(em_data_dir, 'em_168.rds'))
gc()
em_proc(em_filepath = em_266_path, output_filepath = file.path(em_data_dir, 'em_266.rds'))
gc()




# Load DMR calling script
source('src/preprocess/3_make_nano_DMR.R')


# Save phased methylation data for sample '500' for the DMR figure.
h1_500 <- loadHap(file.path(sample_500_dir, get_bedmethyl_files(sample_500_dir)["h1"]))
saveRDS(h1_500, file = file.path(dmrs_out_dir, 's500_h1.rds'))

h2_500 <- loadHap(file.path(sample_500_dir, get_bedmethyl_files(sample_500_dir)["h2"]))
saveRDS(h2_500, file = file.path(dmrs_out_dir, 's500_h2.rds'))

# Run DMR pipeline for multiple samples
dmrs_pipeline(sample_input_dir = sample_51_dir, out_dir = dmrs_out_dir, samplename = "1", min_cov = 5)
gc()

dmrs_pipeline(sample_input_dir = sample_69_dir, out_dir = dmrs_out_dir, samplename = "2", min_cov = 5)
gc()

dmrs_pipeline(sample_input_dir = sample_168_dir, out_dir = dmrs_out_dir, samplename = "3", min_cov = 5)
gc()

dmrs_pipeline(sample_input_dir = sample_266_dir, out_dir = dmrs_out_dir, samplename = "4", min_cov = 5)
gc()


dmrs_pipeline(sample_input_dir = sample_541_dir, out_dir = dmrs_out_dir, samplename = "A", min_cov = 5)
gc()

dmrs_pipeline(sample_input_dir = sample_500_dir, out_dir = dmrs_out_dir, samplename = "B", min_cov = 5)
gc()

dmrs_pipeline(sample_input_dir = sample_577_dir, out_dir = dmrs_out_dir, samplename = "C", min_cov = 5)
gc()

dmrs_pipeline(sample_input_dir = sample_475_dir, out_dir = dmrs_out_dir, samplename = "D", min_cov = 5)
gc()


# Creating Annotation files
source('src/preprocess/4_create_annotations.R')
make_annotations(reference_directory = ref_dir)

