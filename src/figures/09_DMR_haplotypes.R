# ==============================================================
# 09_DMR_haplotypes.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script detects differentially methylated regions (DMRs) between
# two haplotypes using Nanopore sequencing data. It processes methylation
# counts, applies BSSeq statistical tests, and extracts a full DMR.
#
# Dependencies:
# - Requires the DSS and bsseq packages for DMR analysis.
# - Requires methylation data stored in .rds files.
#
# Outputs:
# - Processed BSSeq object.
# - List of detected DMRs.
# - Extracted single full DMR for downstream visualization.
#
# ==============================================================


# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Identify differentially methylated regions (DMRs) between haplotypes
#'
#' This function loads Nanopore methylation data for two haplotypes, filters
#' for CpG sites with sufficient coverage, and applies the BSSeq statistical
#' framework to detect DMRs.
#'
#' @return A data frame containing detected DMRs with genomic coordinates.
dmr_function <- function() {
   set.seed(123)
  # Read sample 500 h1, h2, (sample b)
  nano_500_h1 <- readRDS(file = "data/dmrs/s500_h1.rds")
  nano_500_h2 <- readRDS(file = "data/dmrs/s500_h2.rds")


  min_cov <- 5
  # min_cov <- 10

  nano_500_h1 <- nano_500_h1[nano_500_h1$N >= min_cov, ]
  nano_500_h2 <- nano_500_h2[nano_500_h2$N >= min_cov, ]



  # Make BSSeq data
  BS_500 <- makeBSseqData(
    list(
      nano_500_h1,
      nano_500_h2
    ),
    c("hap1", "hap2")
  )

  BS_500 <- orderBSseq(BS_500, seqOrder = c(paste0("chr", 1:22), "chrX", "chrY"))

  dmlTest.sm.500 <- DMLtest(BS_500, group1 = c("hap1"), group2 = c("hap2"), equal.disp = FALSE, smoothing = TRUE, smoothing.span = 500, ncores = 4)

  dmrs500 <- callDMR(dmlTest.sm.500) # Default p-value threshold is 1e-5


  # Pull a single full DMR


  full <- which.min(abs(dmrs500$start - 53553663))
  full_dmr <- dmrs500[full, ]

  # nano_500_h1
  h1_dmr <- nano_500_h1[(full_dmr$start <= nano_500_h1$pos) &
    (nano_500_h1$pos < full_dmr$end) &
    (full_dmr$chr == nano_500_h1$chr), ]



  h2_dmr <- nano_500_h2[(full_dmr$start <= nano_500_h2$pos) &
    (nano_500_h2$pos < full_dmr$end) &
    (full_dmr$chr == nano_500_h2$chr), ]


  h1_dmr <- h1_dmr %>% mutate(haplotype = "Haplotype 1")
  h2_dmr <- h2_dmr %>% mutate(haplotype = "Haplotype 2")

  full_dmr <- rbind(h1_dmr, h2_dmr)




  # Pull the partial DMR

  partial <- which.min(abs(dmrs500$start - 170243619))
  partial_dmr <- dmrs500[partial, ]


  h1_dmr <- nano_500_h1[(partial_dmr$start <= nano_500_h1$pos) &
    (nano_500_h1$pos < partial_dmr$end) &
    (partial_dmr$chr == nano_500_h1$chr), ]



  h2_dmr <- nano_500_h2[(partial_dmr$start <= nano_500_h2$pos) &
    (nano_500_h2$pos < partial_dmr$end) &
    (partial_dmr$chr == nano_500_h2$chr), ]

  h1_dmr <- h1_dmr %>% mutate(haplotype = "Haplotype 1")
  h2_dmr <- h2_dmr %>% mutate(haplotype = "Haplotype 2")

  partial_dmr <- rbind(h1_dmr, h2_dmr)


  # Make the figures

  dmr_plot_full <- ggplot(data = full_dmr, aes(x = pos, y = (X / N), color = haplotype)) +
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("Haplotype 1" = "goldenrod", "Haplotype 2" = "skyblue")) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor.x = element_blank(),
      #    axis.text.x = element_text(angle = 45, hjust = 1)
      axis.text.x = element_blank()
    ) +
    labs(color = "", x = "", y = "\u03B2 Value")

  # print(dmr_plot_full)


  dmr_plot_partial <- ggplot(data = partial_dmr, aes(x = pos, y = (X / N), color = haplotype)) +
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("Haplotype 1" = "goldenrod", "Haplotype 2" = "skyblue")) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor.x = element_blank(),
      #    axis.text.x = element_text(angle = 45, hjust = 1)
      axis.text.x = element_blank()
    ) +
    labs(color = "", x = "", y = "\u03B2 Value")

  # print(dmr_plot_partial)

  dmr_plot <- dmr_plot_full + dmr_plot_partial + plot_layout(guides = "collect", axis_titles = "collect")
  dmr_plot <- dmr_plot & theme(legend.position = "bottom")


  return(dmr_plot)
}

dmr_plot <- dmr_function()
