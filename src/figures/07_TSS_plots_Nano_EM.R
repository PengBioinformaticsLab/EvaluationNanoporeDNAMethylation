# ==============================================================
# 07_TSS_plots_Nano_EM.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script analyzes CpG sites in relation to transcription start sites (TSS).
# It processes genomic data from Nanopore and EM-Seq datasets, calculates
# distances from CpG sites to TSS, and bins these distances for visualization.
#
# Dependencies:
# - Requires GenomicRanges, IRanges, and tidyverse.
# - Requires reference annotation files (.rds).
#
# Outputs:
# - TSS figures
# ==============================================================

# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Convert a data frame to a GenomicRanges object
#'
#' This function converts a data frame containing genomic coordinates 
#' into a GenomicRanges object for use in TSS distance calculations.
#'
#' @param df Data frame containing columns: chr, start, score, and methylation data.
#' @param nano Binary flag (1 for Nanopore, 0 for EM-Seq) to adjust methylation calculation.
#' @return GRanges object with CpG positions, coverage, and methylation counts
## Convert the dataframes from the .rds files to genomic ranges objects
gr_conv <- function(df, nano = 1) {
  gr_obj <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(
      start = df$start
    ),
    strand = "*", # switch '.' -> '*' for GenomicRanges package.
    coverage = df$score,
    methy = if (nano == 1){
      df$NMod + df$NOtherMod
    } else{
      methy = df$methy
    }
  )
}


#' Calculate distances from CpG sites to the nearest TSS
#'
#' This function computes the distance of CpG sites to the nearest TSS
#' and assigns negative values for upstream sites.
#'
#' @param gr_obj GRanges object containing CpG sites.
#' @param TSS_ref GRanges object containing TSS locations (default: `tss`).
#' @return GRanges object with CpG-TSS distances included in metadata columns.
TSS_carve <- function(gr_obj, TSS_ref = tss) {
  # Find the distance of the nearest TSS for each CpG.
  TSS_dist <- distanceToNearest(gr_obj, TSS_ref)
  TSS_dist <- mcols(TSS_dist)$distance
  
  # Make a 2kb window
  TSS_window <- TSS_dist <= 2000
  TSS_dist <- TSS_dist[TSS_window]
  
  # Find the upstream indices to filter for them within the nearestDistance object.
  TSS_upstream <- precede(gr_obj[TSS_window, ], subject = TSS_ref)
  TSS_near <- nearest(gr_obj[TSS_window, ], TSS_ref)
  upstream_indices <- TSS_upstream == TSS_near
  
  # Multiply upstream indices within the distance object by -1.
  TSS_dist[upstream_indices] <- TSS_dist[upstream_indices] * -1
  
  gr_obj <- gr_obj[TSS_window, ]
  mcols(gr_obj)$distance <- TSS_dist
  
  
  return(gr_obj)
}


#' Bin CpG sites based on distance from TSS
#'
#' This function groups CpG sites into bins based on their proximity
#' to TSS for easier visualization.
#'
#' @param carved_TSS GRanges object with CpG-TSS distances.
#' @return Data frame with CpG counts per distance bin.
TSS_binner <- function(carved_TSS){
  return(
    as.data.frame(carved_TSS) %>%
      dplyr::select(distance, coverage, methy) %>%
      mutate(beta = methy / coverage) %>%
      mutate(bin = cut(distance, breaks = 80)) %>% # (2000bp * 2) / x breaks = 50 bp bins => x = 80
      group_by(bin) %>%
      summarise(
        mean_cov = mean(coverage),
        mean_beta = mean(beta),
        bin_center = mean(as.numeric(sub("\\((.+),(.+)\\]", "\\1", bin)) +
                            as.numeric(sub("\\((.+),(.+)\\]", "\\2", bin))) / 2
      )
  )
}


#Input: binned TSS objects for nanopore and EM-SEQ
TSS_cov_plotter <- function(TSS_nano, TSS_em) {
  
  nanopore_color <- "dodgerblue"
  em_color <- "darkorange"
  
  TSS_nano <- TSS_binner(TSS_nano)
  TSS_em <- TSS_binner(TSS_em)
  
  
  
  TSS_cov_plot <- ggplot() +
    geom_line(data = TSS_nano, aes(x = bin_center, y = mean_cov, color = "Nanopore")) +
    geom_line(data = TSS_em, aes(x = bin_center, y = mean_cov, color = "EM-Seq")) +
    labs(x = "", y = "Mean Coverage", color = "") +
    scale_x_continuous(labels = function(x) ifelse(x == 0, "TSS", paste0(x, "bp"))) +
    scale_y_continuous(limits = c(1, 20), labels = function(y) paste0(y, "X")) +
    scale_color_manual(values = c("EM-Seq" = em_color, "Nanopore" = nanopore_color)) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    )
  
  return(TSS_cov_plot)
}

TSS_beta_plotter <- function(TSS_nano, TSS_em) {
  
  nanopore_color <- "dodgerblue"
  em_color <- "darkorange"
  
  # Bin statistics for each platform at 1X and 10X
  TSS_nano <- TSS_binner(TSS_nano)
  TSS_em <- TSS_binner(TSS_em)
  
  
  TSS_cov_plot <- ggplot() + 
    geom_line(data = TSS_nano,  aes(x = bin_center, y = mean_beta, color = "Nanopore")) +
    geom_line(data = TSS_em, aes(x = bin_center, y = mean_beta, color = "EM-Seq")) +
    scale_x_continuous(labels = function(x) ifelse(x == 0, "TSS", paste0(x, "bp"))) +
    scale_y_continuous(limits = c(0, 1), labels = function(y) round(y, 2)) +
    labs(x ="" , y = "Methylation Ratio", color = "" ) + 
    scale_color_manual(values = c("EM-Seq" = em_color, "Nanopore" = nanopore_color)) +
    theme_minimal()  + 
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    )
  #    theme(legend.position = "none")
  
  return(TSS_cov_plot)
}




count_ratio_plotter <- function(TSS_nano, TSS_em) {
  
  #### 1X coverage minimum #######
  # Count CpG sites, grouped by distance from TSS for each method:
  nano_count_1X <- as.data.frame(TSS_nano) %>%
    group_by(distance) %>%
    summarise(nano_cpg_count = n())
  
  em_count_1X <- as.data.frame(TSS_em) %>%
    group_by(distance) %>%
    summarise(em_cpg_count = n())
  
  
  # Join the two TSS objects together by distance
  counts_by_distance_1X <- inner_join(
    x = nano_count_1X,
    y = em_count_1X,
    by = "distance",
    suffix = c("nano", "em")
  )
  
  
  #### 10X coverage minimum #######
  nano_count_10X <- as.data.frame(TSS_nano) %>%
    filter(coverage >= 10) %>%
    group_by(distance) %>%
    summarise(nano_cpg_count_10X = n())
  
  em_count_10X <- as.data.frame(TSS_em) %>%
    filter(coverage >= 10) %>%
    group_by(distance) %>%
    summarise(em_cpg_count_10X = n())
  
  counts_by_distance_10X <- inner_join(
    x = nano_count_10X,
    y = em_count_10X,
    by = "distance",
    suffix = c("nano", "em")
  )
  
  #### Brining it all together #### 
  counts_by_distance <- inner_join(
    x = counts_by_distance_1X,
    y = counts_by_distance_10X,
    by = "distance"
  )
  
  counts_by_distance <- counts_by_distance %>%
    mutate(
      count_ratio_1X = em_cpg_count / nano_cpg_count,
      count_ratio_10X = em_cpg_count_10X / nano_cpg_count_10X
    )
  
  
  ### Plot
  
  color_1X <- "salmon"
  color_10X <- "skyblue"
  
  ratio_2kb <- ggplot(data = counts_by_distance) +
    #geom_line(aes(x = distance, y = count_ratio_1X, color = "1X Coverage")) +
    #geom_line(aes(x = distance, y = count_ratio_10X, color = "10X Coverage")) +
    geom_line(aes(x = distance, y = count_ratio_1X, color = "≥ 1x")) +
    geom_line(aes(x = distance, y = count_ratio_10X, color = "≥ 10x")) +
    labs(x = "", y = "CpG Count Ratio \n (EM-Seq / Nanopore)", color = "", title = "") +
    scale_x_continuous(labels = function(x) ifelse(x == 0, "TSS", paste0(x, "bp"))) +
    scale_y_continuous(limits = c(0, 1), labels = function(x) sprintf("%.2f", x)) +
    scale_color_manual(values = c("≥ 1x" = color_1X, "≥ 10x" = color_10X)) +
    #scale_color_manual(values = c("1X Coverage" = color_1X, "10X Coverage" = color_10X)) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    )
  
  
  
  return(ratio_2kb)
}


hype_plotter <- function(TSS_df) {
  #### df-Seq hypo/hyper lines
  
  dichom_df <- as.data.frame(TSS_df) %>%
    mutate(beta = methy * (100 / coverage)) %>%
    filter(beta <= 10 | beta >= 90) %>%
    mutate(hypo = (beta <= 10))
  
  
  df_lines_hypo <- dichom_df %>%
    filter(hypo == 1) %>%
    mutate(bin = cut(distance, breaks = 80)) %>% # 50-bp bins
    group_by(bin) %>%
    summarise(
      mean_cov = mean(coverage),
      mean_beta = mean(beta),
      bin_center = mean(as.numeric(sub("\\((.+),(.+)\\]", "\\1", bin)) +
                          as.numeric(sub("\\((.+),(.+)\\]", "\\2", bin))) / 2
    )
  
  df_lines_hyper <- dichom_df %>%
    filter(hypo == 0) %>%
    mutate(bin = cut(distance, breaks = 80)) %>% # 50-bp bins
    group_by(bin) %>%
    summarise(
      mean_cov = mean(coverage),
      mean_beta = mean(beta),
      bin_center = mean(as.numeric(sub("\\((.+),(.+)\\]", "\\1", bin)) +
                          as.numeric(sub("\\((.+),(.+)\\]", "\\2", bin))) / 2
    )
  
  ### Plotting
  
  hyper_color <-"forestgreen" 
  hypo_color <- "#957DAD"
  
  
  hype_plot <- ggplot() +
    geom_line(data = df_lines_hyper, aes(x = bin_center, y = mean_cov, color = "High Methylation (\u03B2 \u2265 0.9)")) +
    geom_line(data = df_lines_hypo, aes(x = bin_center, y = mean_cov, color = "Low Methylation (\u03B2 \u2264 0.1)")) +
    labs(x = "", y = "Mean Coverage", color = "", title = "") +
    scale_x_continuous(labels = function(x) ifelse(x == 0, "TSS", paste0(x, "bp"))) +
    scale_y_continuous(limits = c(1, 20), labels = function(y) paste0(y, "X")) +
    scale_color_manual(values = c("High Methylation (\u03B2 \u2265 0.9)" = hyper_color, "Low Methylation (\u03B2 \u2264 0.1)" = hypo_color)) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom")
  
  hype_plot <- hype_plot
  
  
  
  return(hype_plot)
}




TSS_beta_plotter_new <- function(TSS_1X, TSS_10X, color = "nano") {
  # Define colors and labels based on platform
  color_palette <- list(
    nano = list(
      color_1X = "dodgerblue", color_10X = "dodgerblue3",
      all_cpgs = "Nanopore (≥ 1x)", atleast10 = "Nanopore (≥ 10x)"
    ),
    emseq = list(
      color_1X = "darkorange", color_10X = "darkorange3",
      all_cpgs = "EM-Seq (≥ 1x)", atleast10 = "EM-Seq (≥ 10x)"
    )
  )
  
  TSS_1X <- TSS_binner(TSS_1X)
  TSS_10X <- TSS_binner(TSS_10X)

  # Select color scheme based on input
  selected_colors <- if (color == "nano") color_palette$nano else color_palette$emseq

  # Generate the plot
  TSS_cov_plot <- ggplot() +
    geom_line(data = TSS_1X, aes(x = bin_center, y = mean_beta, color = selected_colors$all_cpgs), linetype = 1) +
    geom_line(data = TSS_10X, aes(x = bin_center, y = mean_beta, color = selected_colors$atleast10), linetype = 2) +
    scale_x_continuous(labels = function(x) ifelse(x == 0, "TSS", paste0(x, "bp"))) +
    scale_y_continuous(limits = c(0, 1), labels = function(x) sprintf("%.2f", x)) +
#    labs(x = "", y = "Methylation\nProportion", color = "") +
    labs(x = "", y = "Mean \u03B2 Value", color = "") +
    scale_color_manual(values = setNames(
      c(selected_colors$color_1X, selected_colors$color_10X),
      c(selected_colors$all_cpgs, selected_colors$atleast10)
    )) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.text = element_text(size = 8)
    )

  return(TSS_cov_plot)
}


# ==============================================================
# EXECUTE ANALYSIS
# ==============================================================

#### Load TSS reference:


gene_annotations <- readRDS(file = 'data/reference/gene_annotations.rds')
# pull the transcripts only
tss <- gene_annotations[gene_annotations$type == "transcript"]
# Shrink down to just the start site.
ranges(tss) <- IRanges(start = start(tss), end = start(tss) + 1)


#### Creating plots



make_TSS_figure <- function(nano_df, em_df) {
  ## data preprocessing
  nano_gr <- gr_conv(nano_df)
  nano_TSS <- TSS_carve(gr_obj = nano_gr, TSS_ref = tss)

  em_gr <- gr_conv(em_df, nano = 0)
  em_TSS <- TSS_carve(gr_obj = em_gr, TSS_ref = tss)

  ## Coverage plot
  coverage_plot <- TSS_cov_plotter(nano_TSS, em_TSS)

  ## Methylation plots
  hype_nano <- hype_plotter(TSS_df = nano_TSS) + labs(title = "Nanopore")
  hype_em <- hype_plotter(TSS_df = em_TSS) + labs(title = "EM-Seq")

  ## Count ratio plot
  count_plot <- count_ratio_plotter(nano_TSS, em_TSS)
  
  TSS_beta_nano <- TSS_beta_plotter_new(TSS_1X = nano_TSS, TSS_10X = nano_TSS[nano_TSS$coverage >= 10, ], color = "nano")
  TSS_beta_em <- TSS_beta_plotter_new(TSS_1X = em_TSS, TSS_10X = em_TSS[em_TSS$coverage >= 10, ], color = "em")
  
  beta_p <- TSS_beta_nano +
    TSS_beta_em +
    plot_layout(axis_title = "collect") +
    plot_layout(axis_title = "collect", guides = "collect") & theme(legend.position = "bottom")
  

  ## All together now:

  layout <- "
AB
CC
DD
"

  hype <- hype_nano + hype_em +
    plot_layout(axis_title = "collect", guides = "collect") & theme(legend.position = "bottom")


  tasteful_TSS <- coverage_plot + count_plot + hype + beta_p + 
    plot_layout(design = layout, axis_titles = "collect") & theme(plot.margin = margin(14, 10, 14, 10))

  tasteful_TSS <- tasteful_TSS + plot_annotation(tag_levels = "A")

  return(tasteful_TSS)
}







ggsave(filename = 'output/TSS_sample1.png', plot=make_TSS_figure(nano_1, em_1), width = 8, height = 10)
ggsave(filename = 'output/TSS_sample2.png', plot=make_TSS_figure(nano_2, em_2), width = 8, height = 10)
ggsave(filename = 'output/TSS_sample3.png', plot=make_TSS_figure(nano_3, em_3), width = 8, height = 10)
ggsave(filename = 'output/TSS_sample4.png', plot=make_TSS_figure(nano_4, em_4), width = 8, height = 10)

