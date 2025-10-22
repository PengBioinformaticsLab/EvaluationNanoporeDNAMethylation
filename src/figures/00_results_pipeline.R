# ==============================================================
# 00_results_pipeline.R
# Author: Steven Brooks
# Date: 02/14/25
# --------------------------------------------------------------
# Description:
# This script generates various figures comparing Nanopore and 
# EM-Seq methylation data. It sources external scripts for heatmaps,
# dotplots, violin plots, TSS plots, and DMR statistics. It also
# loads required datasets and saves final figures.
#
# Dependencies:
# - Requires ggplot2, tidyverse, ggpubr, viridis, GenomicRanges, etc.
# - Loads dataset paths from '01_filepaths.R'
# - Calls multiple figure generation scripts
#
# Outputs:
# - Saves figures in 'output/' directory
#
# ===============================================================

# Set seed for reproducibility
set.seed(123)

# Load necessary libraries for visualization and statistical analysis
library(ggplot2)    # Data visualization
library(tidyverse)  # Data manipulation
library(ggpubr)     # Publication-ready plots
library(viridis)    # Color scales
library(GenomicRanges) # Genomic range handling
library(gghalves)   # Half violin plots
library(patchwork)  # Plot arrangement
library(DSS)        # Differential methylation statistics
library(bsseq)      # Bisulfite sequencing analysis
library(testthat)   # Code validation



#' Safely Source an R Script
#'
#' This function checks whether a specified R script exists before sourcing it. 
#' If the script is missing, it issues a warning instead of throwing an error.
#'
#' @param script_path Character string specifying the path to the script.
#' @return Invisibly returns `NULL`. If the script exists, it is sourced; otherwise, a warning is issued.
source_script <- function(script_path) {
  if (file.exists(script_path)) {
    source(script_path)
  } else {
    warning("Missing script: ", script_path)
  }
}


# ==============================================================
# DATA LOADING
# ==============================================================

# Load file paths for datasets (Defined in 01_filepaths.R)
source_script('src/figures/01_filepaths.R')


# ==============================================================
# FIGURE GENERATION
# ==============================================================

# Generate heatmaps for EM-Seq at different coverage levels
source_script('src/figures/02_heatmap.R')

# Generate summary table comparing Nanopore and EM-Seq
source_script('src/figures/03_stats_Nano_EM.R')

# Generate summary table for 5mC and 5hmC modtypes
source_script('src/figures/04_modtype_stats.R')

# Generate dotplots for Nanopore vs. EM-Seq methylation levels
source_script('src/figures/05_dotplot_Nano_EM.R')

# Generate violin plots to compare distributions of methylation coverage
source_script('src/figures/06_violin_plot_Nano_EM.R')



# ==============================================================
# FIGURE 2: DOTPLOT + VIOLIN PLOT
# ==============================================================

### First RED memory error, reran 01 an start from here:

# Generate a combined plot with dotplot and violin plot
# Combines individual plot functions from source_scriptd scripts
figure_2 <- (make_dotplot(cpg_table) + theme(legend.position= "bottom")) + 
  violin_plotter(nano_1, em_1) + plot_layout(guides = 'collect') & 
  theme(legend.justification = "center", legend.position = "bottom") 

figure_2 <- figure_2 + plot_annotation(tag_levels = 'A')

ggsave(filename = "output/figure2_dotplot_violins.png", plot = figure_2, dpi = 300, width = 16)

# ==============================================================
# TSS PLOTS: NANOPORE VS EM-SEQ
# ==============================================================

##### TSS plots (Nanopore vs. EM-Seq)
source_script('src/figures/07_TSS_plots_Nano_EM.R')

# ==============================================================
# MODIFICATION TYPE ANALYSIS
# ==============================================================

# Load script for modification type density plots
source_script('src/figures/08_modtype_densities.R')


# ==============================================================
# DMR  + Figure 4
# ==============================================================


### Red 

##### DMR between haplotypes 
source_script('src/figures/09_DMR_haplotypes.R')

#Generate Figure 4:  Fig A,B modtype for blood and brain / Fig C,D  full and partial DMRs
methyl_plot <- methyl_dplot(nano_a) + methyl_dplot(nano_1) + plot_layout(guides = 'collect', axis_titles = 'collect') & theme(legend.position = "bottom")

figure4 <- methyl_plot / dmr_plot + plot_layout(heights = c(4,4))
figure4 <- figure4 + plot_annotation(tag_levels = 'A')

ggsave('output/figure4.png', plot = figure4, dpi = 300, width = 12)

# RED

##### DMR Statistics
source_script('src/figures/10_DMR_stats.R')

