# ==============================================================
# 05_dotplot_Nano_EM.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script counts CpG sites at different coverage thresholds
# for both Nanopore and EM-Seq methylation datasets. It then 
# formats the results into a table for visualization in a dot plot.
#
# Dependencies:
# - Requires input data frames from Nanopore and EM-Seq (.rds files).
# - Requires ggplot2 for visualization.
#
# Outputs:
# - Returns a table of CpG site counts at different coverage levels.
#
# ==============================================================


# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Count CpG sites at different coverage thresholds
#'
#' This function iterates through coverage thresholds (1X to 30X)
#' and counts the number of CpG sites exceeding each threshold.
#'
#' @param df Data frame containing CpG methylation data.
#' @return Numeric vector of CpG counts at each coverage level.
cpg_counter <- function(df){
  
  if (!"score" %in% colnames(df)) {
    stop("Column 'score' not found in the input data frame")
  }
  
  # Define maximum coverage threshold
  max_coverage <- 30
  
  # Initialize vector to store CpG counts
  numCpG_Cov <- numeric(max_coverage)
  
  # Iterate through coverage levels and count CpGs above threshold
  for (i in 1:max_coverage){
    numCpG_Cov[i] <- sum(df$score >= i)
  }
  
  return(numCpG_Cov)
}


#' Create a table of CpG counts for Nanopore and EM-Seq datasets
#'
#' This function constructs a data frame containing the number of
#' detected CpG sites at each coverage threshold for multiple samples.
#'
#' @param nano_1, nano_2, nano_3, nano_4 Data frames for Nanopore samples.
#' @param em_1, em_2, em_3, em_4 Data frames for EM-Seq samples.
#' @return Data frame with CpG counts across coverage levels.
create_cpg_table <- function(nano_1, em_1,
                             nano_2, em_2,
                             nano_3, em_3,
                             nano_4, em_4) {
  
  # Initialize coverage levels (1X to 30X)
  cpg_count <- c(1:30)
  
  # Compute CpG counts for each dataset
  cpg_count <- cbind(cpg_count, cpg_counter(nano_1))
  cpg_count <- cbind(cpg_count, cpg_counter(nano_2))
  cpg_count <- cbind(cpg_count, cpg_counter(nano_3))
  cpg_count <- cbind(cpg_count, cpg_counter(nano_4))
  cpg_count <- cbind(cpg_count, cpg_counter(em_1))
  cpg_count <- cbind(cpg_count, cpg_counter(em_2))
  cpg_count <- cbind(cpg_count, cpg_counter(em_3))
  cpg_count <- cbind(cpg_count, cpg_counter(em_4))
  
  # Convert to data frame format
  cpg_count <- data.frame(cpg_count)
  
  # Assign meaningful column names
  colnames(cpg_count) <- c(
    "coverage",
    "Nanopore S1",
    "Nanopore S2",
    "Nanopore S3",
    "Nanopore S4",
    "EM S1",
    "EM S2",
    "EM S3",
    "EM S4"
  )
  
  return(cpg_count)
}


# ==============================================================
# GENERATE DOT PLOT DATA
# ==============================================================

#' Generate a dot plot visualizing CpG coverage distributions
#'
#' This function generates a dot plot comparing CpG site coverage
#' distributions across multiple samples for Nanopore and EM-Seq.
#'
#' @param cpg_table Data frame containing CpG counts across coverage levels.
#' @return ggplot2 object representing the dot plot.
make_dotplot <- function(cpg_table) {
  
  base_font_size = 12
  
  colors <- c(
    "Nanopore S1" = "lightblue", "EM S1" = "#FFD580",
    "Nanopore S2" = "skyblue", "EM S2" = "gold",
    "Nanopore S3" = "dodgerblue", "EM S3" = "darkorange",
    "Nanopore S4" = "#4169E1", "EM S4" = "goldenrod",
    "Nanopore" = "dodgerblue", "EM-Seq" = "darkorange"
  )
  
  dummy_data_legend <- data.frame(coverage = c(-1,-1),
                                  y = c(-1, -1))
  
  cpg_samples <- ggplot() +
    geom_point(data = cpg_table, aes(x = coverage, y = `EM S1`, color = "EM S1"), show.legend = FALSE) +
    geom_point(data = cpg_table, aes(x = coverage, y = `Nanopore S1`, color = "Nanopore S1"), show.legend = FALSE) +
    geom_point(data = cpg_table, aes(x = coverage, y = `EM S2`, color = "EM S2"), show.legend = FALSE) +
    geom_point(data = cpg_table, aes(x = coverage, y = `Nanopore S2`, color = "Nanopore S2"), show.legend = FALSE) +
    geom_point(data = cpg_table, aes(x = coverage, y = `EM S3`, color = "EM S3"), show.legend = FALSE) +
    geom_point(data = cpg_table, aes(x = coverage, y = `Nanopore S3`, color = "Nanopore S3"), show.legend = FALSE) +
    geom_point(data = cpg_table, aes(x = coverage, y = `EM S4`, color = "EM S4"), show.legend = FALSE) +
    geom_point(data = cpg_table, aes(x = coverage, y = `Nanopore S4`, color = "Nanopore S4"), show.legend = FALSE) +
    geom_point(data = dummy_data_legend, aes(x = coverage, y = y, color = "Nanopore")) + 
    geom_point(data = dummy_data_legend, aes(x = coverage, y = y, color = "EM-Seq")) + 
    # scale_color_manual(
    #   values = colors,
    #    breaks = c(
    #      "Nanopore S1", "Nanopore S2", "Nanopore S3", "Nanopore S4", 
    #      "EM S1","EM S2", "EM S3", "EM S4"))+
    scale_color_manual(
      name = "",
      # Show only these two levels in the legend:
      breaks = c("Nanopore", "EM-Seq"),
      # Rename them in the legend:
      labels = c("Nanopore" = "Nanopore", "EM-Seq" = "EM‑Seq"),
      # Assign actual plot colors to both Sample and Platform levels
      values = colors
    ) + 
    labs(y = "Number of CpG Sites", x = "Minimum Coverage Threshold", color = "") +
    scale_y_continuous(
      limits = c(1e6, 30e6), 
      breaks = c(1e6, 15e6, 30e6),
      labels = function(y) paste0(y / 1e6, "M")
    ) +
    scale_x_continuous(
      limits = c(1,30),
      breaks = c(1, 10, 20, 30),
      labels = function(x) paste0(x, "X")
    ) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(size = base_font_size),               # X-axis labels
      axis.text.y = element_text(size = base_font_size),               # Y-axis labels
      axis.title.x = element_text(size = base_font_size + 2),          # X-axis title (+2 adjustment)
      axis.title.y = element_text(size = base_font_size + 2),          # Y-axis title (+2 adjustment)
      legend.text = element_text(size = base_font_size),               # Legend labels
      legend.title = element_text(size = base_font_size + 2)           # Legend title (+2 adjustment)
    )
  
  
  
  
  return(cpg_samples)
}


# ==============================================================
# EXECUTE ANALYSIS
# ==============================================================

# Generate CpG count table
cpg_table <- create_cpg_table(nano_1, em_1, 
                              nano_2, em_2, 
                              nano_3, em_3, 
                              nano_4, em_4)





# Generate dot plot
dotplot <- make_dotplot(cpg_table)

# Save plot to output directory
ggsave("output/dotplot.png", plot = dotplot, dpi = 300, width = 8, height = 6)

# ==============================================================
# UNIT TESTING
# ==============================================================

# Create mock datasets for testing
test_df <- data.frame(
  Coverage = c(5, 10, 15, 20, 25, 30),
  CpG_Count = c(100, 80, 60, 40, 20, 10)
)

empty_df <- data.frame(Coverage = numeric(0), CpG_Count = numeric(0))

# Test cases for cpg_counter

test_that("cpg_counter handles missing columns", {
  expect_error(cpg_counter(data.frame(OtherColumn = c(1, 2, 3))))
})

# Test cases for make_dotplot

test_that("make_dotplot creates a ggplot object", {
  plot <- make_dotplot(test_df)
  expect_true("gg" %in% class(plot))
})






