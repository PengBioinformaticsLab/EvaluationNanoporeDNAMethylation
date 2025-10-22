# ==============================================================
# 06_violin_plot_Nano_EM.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script generates violin plots to visualize the distribution 
# of CpG coverage levels across Nanopore and EM-Seq datasets.
# It also processes genomic ranges and overlaps with gene annotations.
#
# Dependencies:
# - Requires ggplot2, GenomicRanges, and IRanges for visualization.
# - Requires reference annotation files (.rds).
#
# Outputs:
# - Violin plots saved as PNG images in the output directory.
#
# ==============================================================


# ==============================================================
# LOAD REFERENCE DATA
# ==============================================================

# Load CpG annotation (islands, shores, etc.)
cpg_annotation <- readRDS(file = 'data/reference/cpg_annotation.rds')

# Load gene annotations (hg38 reference)
gene_annotations <- readRDS(file = 'data/reference/gene_annotations.rds')


# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Convert a data frame to a GenomicRanges object
#'
#' This function converts a data frame containing genomic coordinates 
#' into a GenomicRanges object for use with the Bioconductor ecosystem.
#'
#' @param df Data frame containing columns: chr, start, score.
#' @return GRanges object with CpG site positions and coverage values.
## Convert the dataframes from the .rds files to genomic ranges objects
gr_conv <- function(df) {
  gr_obj <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(
      start = df$start
    ),
    strand = "*", # switch '.' -> '*' for GenomicRanges package.
    coverage = df$score
  )
}


#' Extract CpG sites from a specific genomic region
#'
#' Given a GenomicRanges object and gene annotations, this function 
#' returns unique CpG sites located within the specified region.
#'
#' @param gr_obj GRanges object containing CpG sites.
#' @param gene_annotations GRanges object with gene regions.
#' @return Data frame containing unique CpG sites within the region.
carve_region <- function(gr_obj, gene_annotations) {

  overlaps <- findOverlaps(gr_obj, gene_annotations)
  overlapping_gr <- gr_obj[queryHits(overlaps)]

  overlapping_gr  <- unique(overlapping_gr)

  return(as.data.frame(overlapping_gr))
}


#' Generate a violin plot of CpG coverage distributions
#'
#' This function creates a violin plot comparing the coverage distributions
#' of Nanopore and EM-Seq methylation datasets.
#'
#' @param nano_df Data frame containing Nanopore CpG methylation data.
#' @param em_df Data frame containing EM-Seq CpG methylation data.
#' @return ggplot2 object representing the violin plot.
violin_plotter <- function(nano_df, em_df) {


  gr_nano <- gr_conv(nano_df)
  gr_em <- gr_conv(em_df)

  # Figure options
  my_nudge <- .01
  my_bw <- 1
  base_font_size <- 12

  my_y_lim <- 30 # round(max(mean(gr_nano$coverage), mean(gr_em$coverage)) * (3/5)) * 5

  #   unique(gene_annotations$type)
  # Levels: gene transcript exon CDS start_codon stop_codon UTR


  # Truncate excessively read CpG sites to a maximum coverage
  max_cov <- 100
  gr_em$coverage <- ifelse(gr_em$coverage > max_cov, max_cov + 1, gr_em$coverage)
  gr_nano$coverage <- ifelse(gr_nano$coverage > max_cov, max_cov + 1, gr_nano$coverage)



  vp <- ggplot() +
    ######### Genes
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = gene_annotations[gene_annotations$type == "gene"]
      ),
      aes(x = factor("Genes"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = gene_annotations[gene_annotations$type == "gene"]
      ),
      aes(x = factor("Genes"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### transcript
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = gene_annotations[gene_annotations$type == "transcript"]
      ),
      aes(x = factor("Transcripts"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = gene_annotations[gene_annotations$type == "transcript"]
      ),
      aes(x = factor("Transcripts"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### exon
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = gene_annotations[gene_annotations$type == "exon"]
      ),
      aes(x = factor("Exons"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = gene_annotations[gene_annotations$type == "exon"]
      ),
      aes(x = factor("Exons"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### CDS
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = gene_annotations[gene_annotations$type == "CDS"]
      ),
      aes(x = factor("CDS"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = gene_annotations[gene_annotations$type == "CDS"]
      ),
      aes(x = factor("CDS"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### UTR
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = gene_annotations[gene_annotations$type == "UTR"]
      ),
      aes(x = factor("UTR"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = gene_annotations[gene_annotations$type == "UTR"]
      ),
      aes(x = factor("UTR"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### CpG_Island
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = cpg_annotation[cpg_annotation$type == "CpG_Island"]
      ),
      aes(x = factor("CpG Island"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = cpg_annotation[cpg_annotation$type == "CpG_Island"]
      ),
      aes(x = factor("CpG Island"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### Shore
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = cpg_annotation[cpg_annotation$type == "Shore"]
      ),
      aes(x = factor("Shore"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = cpg_annotation[cpg_annotation$type == "Shore"]
      ),
      aes(x = factor("Shore"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### Shelf
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = cpg_annotation[cpg_annotation$type == "Shelf"]
      ),
      aes(x = factor("Shelf"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = cpg_annotation[cpg_annotation$type == "Shelf"]
      ),
      aes(x = factor("Shelf"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### Open_sea
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = cpg_annotation[cpg_annotation$type == "Open_Sea" ]
      ),
      aes(x = factor("Open Sea"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = cpg_annotation[cpg_annotation$type == "Open_Sea" ]
      ),
      aes(x = factor("Open Sea"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### Promoter
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_nano,
        gene_annotations = cpg_annotation[cpg_annotation$type == "Promoter"]
      ),
      aes(x = factor("Promoter"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = carve_region(
        gr_obj = gr_em,
        gene_annotations = cpg_annotation[cpg_annotation$type == "Promoter"]
      ),
      aes(x = factor("Promoter"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######## Intergenic
    geom_half_violin(
      data = as.data.frame(unique(gr_nano[-queryHits(findOverlaps(gr_nano, gene_annotations[gene_annotations$type == 'gene']))])),
      aes(x = factor("Intergenic"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = as.data.frame(unique(gr_em[-queryHits(findOverlaps(gr_em, gene_annotations[gene_annotations$type == 'gene']))])),
      aes(x = factor("Intergenic"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### Whole Genome
    geom_half_violin(
      data = as.data.frame(gr_nano),
      aes(x = factor("Whole Genome"), y = coverage, color = "Nanopore"),
      side = "l", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    geom_half_violin(
      data = as.data.frame(gr_em),
      aes(x = factor("Whole Genome"), y = coverage, color = "EM-Seq"),
      side = "r", nudge = my_nudge, bw = my_bw, trim = TRUE, show.legend = FALSE
    ) +
    ######### Figure Options
    labs(
      y = "Coverage",
      x = "Genomic Region",
      color = ""
    ) +
    scale_color_manual(values = c("EM-Seq" = "darkorange", "Nanopore" = "dodgerblue")) +
    scale_y_continuous(
      breaks = c(1, 10, 20, 30),
      labels = function(y) paste0(y, "X")
    ) +
    coord_cartesian(ylim = c(0, my_y_lim)) +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(angle = 45, hjust = .75, size = base_font_size), #changed 12/7/24
 #     axis.text.x = element_text(angle = 45, hjust = .5, size = base_font_size),
      axis.text.y = element_text(size = base_font_size),               # Y-axis labels
      axis.title.x = element_text(size = base_font_size + 2),          # X-axis title (+2 adjustment)
      axis.title.y = element_text(size = base_font_size + 2),          # Y-axis title (+2 adjustment)
    )

  return(vp)
}

# ==============================================================
# PLOT
# ==============================================================


ggsave(filename = "output/violins_sample1.png", plot = violin_plotter(nano_1, em_1), width = 12, height = 5)
ggsave(filename = "output/violins_sample2.png", plot = violin_plotter(nano_2, em_2), width = 12, height = 5)
ggsave(filename = "output/violins_sample3.png", plot = violin_plotter(nano_3, em_3), width = 12, height = 5)
ggsave(filename = "output/violins_sample4.png", plot = violin_plotter(nano_4, em_4), width = 12, height = 5)


