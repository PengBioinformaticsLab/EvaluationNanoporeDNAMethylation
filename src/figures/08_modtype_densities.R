# ==============================================================
# 08_modtype_densities.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script generates density plots for methylation modifications 
# (5mC, 5hmC, and combined modifications) using Nanopore sequencing data.
# The density plots visualize the distribution of methylation beta values.
#
# Dependencies:
# - Requires ggplot2 for visualization.
# - Input should be a bedmethyl-like dataframe with methylation data.
#
# Outputs:
# - Returns a ggplot object containing density curves for methylation types.
#
# ==============================================================


# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Generate density plot of methylation modification types
#'
#' This function creates a density plot for three types of methylation 
#' modifications: 5mC, 5hmC, and their combined values.
#'
#' @param nano_bed Data frame containing Nanopore methylation data.
#'                 Expected columns: NMod (5mC count), NOtherMod (5hmC count), score (coverage).
#' @param my_bw Numeric. Bandwidth parameter for smoothing (default = 0.06).
#' @return ggplot2 object with density curves for 5mC, 5hmC, and combined modifications.
methyl_dplot <- function(nano_bed, my_bw = 0.06) {
  
  # Define base font size for plot labels
  base_font_size <- 12
  
  # Generate density plot
  dens_plt <- ggplot(data = nano_bed) +
    
    # Combined 5mC + 5hmC density curve
    stat_density(
      aes(x = (NMod + NOtherMod) / score, color = "5mC and 5hmC"),
      geom = "line",
      position = "identity",
      bw = my_bw
    ) +
    
    # 5mC density curve
    stat_density(
      aes(x = NMod / score, color = "5mC"),
      geom = "line",
      position = "identity",
      bw = my_bw
    ) +
    
    # 5hmC density curve
    stat_density(
      aes(x = NOtherMod / score, color = "5hmC"),
      geom = "line",
      position = "identity",
      bw = my_bw
    ) +
    
    # Define colors for each density curve
    scale_color_manual(values = c("5mC and 5hmC" = "#000080", 
                                  "5mC" = "#228B22", 
                                  "5hmC" = "#E0115F")) +
    
    # Define axis labels and legend
    labs(x = "\u03B2 Value", y = "Density", title = "", color = "") +
    
    # Adjust legend position
    theme(legend.position = "bottom") +
    
    # Use minimal theme with white background
    theme_minimal() +
    
    # Customize plot appearance
    theme(
      panel.border = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Customize text sizes
      axis.text.x = element_text(size = base_font_size), # X-axis labels
      axis.text.y = element_text(size = base_font_size), # Y-axis labels
      axis.title.x = element_text(size = base_font_size + 2), # X-axis title (+2 adjustment)
      axis.title.y = element_text(size = base_font_size + 2), # Y-axis title (+2 adjustment)
      legend.text = element_text(size = base_font_size), # Legend labels
      legend.title = element_text(size = base_font_size + 2) # Legend title (+2 adjustment)
    )
  
  return(dens_plt)
}


# ==============================================================
# EXECUTE DENSITY PLOT GENERATION (EXAMPLE USAGE)
# ==============================================================

# Example usage (uncomment to generate a plot)
# d <- methyl_dplot(nano_1)
#ggsave("output/methyl_density_plot.png", plot = d, dpi = 300, width = 8, height = 6)



## Legacy code below:
# methyl_dplot <- function(nano_bed, my_bw = .06) {
#   # input: 1 nano .bedmethyl dataframe
#   # optionally: a smoothing value
#   # output: density curve figure
#   base_font_size <- 12
# 
#   dens_plt <- ggplot(data = nano_bed) +
#     stat_density(
#       aes(x = (NMod + NOtherMod) / score, color = "5mC and 5hmC"),
#       geom = "line",
#       position = "identity",
#       bw = my_bw
#     ) +
#     stat_density(
#       aes(x = NMod / score, color = "5mC"),
#       geom = "line",
#       position = "identity",
#       bw = my_bw
#     ) +
#     stat_density(
#       aes(x = NOtherMod / score, color = "5hmC"),
#       geom = "line",
#       position = "identity",
#       bw = my_bw
#     ) +
#     scale_color_manual(values = c("5mC and 5hmC" = "#000080", "5mC" = "#228B22", "5hmC" = "#E0115F")) +
#     labs(x = "\u03B2 Value", y = "Density", title = "", color = "") +
#     theme(legend.position = "bottom") +
#     theme_minimal() +
#     theme(
#       panel.border = element_blank(),
#       panel.grid.minor.x = element_blank(),
#       plot.background = element_rect(fill = "white", color = NA),
#       panel.background = element_rect(fill = "white", color = NA),
#       axis.text.x = element_text(size = base_font_size), # X-axis labels
#       axis.text.y = element_text(size = base_font_size), # Y-axis labels
#       axis.title.x = element_text(size = base_font_size + 2), # X-axis title (+2 adjustment)
#       axis.title.y = element_text(size = base_font_size + 2), # Y-axis title (+2 adjustment)
#       legend.text = element_text(size = base_font_size), # Legend labels
#       legend.title = element_text(size = base_font_size + 2) # Legend title (+2 adjustment)
#     )
# 
# 
#   return(dens_plt)
# }



#d <- methyl_dplot(nano_1)

# ggsave(filename = "output/modtype_curve_sample1.png", plot = methyl_dplot(nano_1), width = 8, height = 5)
# ggsave(filename = "output/modtype_curve_sample2.png", plot = methyl_dplot(nano_2), width = 8, height = 5)
# ggsave(filename = "output/modtype_curve_sample3.png", plot = methyl_dplot(nano_3), width = 8, height = 5)
# ggsave(filename = "output/modtype_curve_sample4.png", plot = methyl_dplot(nano_4), width = 8, height = 5)
