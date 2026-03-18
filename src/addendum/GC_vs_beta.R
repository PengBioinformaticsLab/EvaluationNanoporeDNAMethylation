# Load libraries
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(stringr)

out_dir <- "addendum/res"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ONT files
file_paths_ONT <- list.files(path = "data/nanopore", pattern = "nano_.*\\.rds$", full.names = TRUE)
file_order_ONT <- order(as.numeric(str_extract(basename(file_paths_ONT), "\\d+")))
file_paths_ONT <- file_paths_ONT[file_order_ONT]

# EM files
file_paths_EM <- list.files(path = "data/em-seq", pattern = "em_.*\\.rds$", full.names = TRUE)
file_order_EM <- order(as.numeric(str_extract(basename(file_paths_EM), "\\d+")))
file_paths_EM <- file_paths_EM[file_order_EM]

# 1. Load the Reference Tile Data
message("Loading 1kb tile annotations...")
nuc_dt <- fread("addendum/data/tiles_1kb_nuc.bed", 
                select = c(1, 2, 3, 5), 
                col.names = c("chr", "start", "end", "pct_gc"))

# Add a row ID to act as our tile_id identifier
nuc_dt[, tile_id := .I]

gr_tiles <- GRanges(seqnames = nuc_dt$chr,
                    ranges = IRanges(start = nuc_dt$start, end = nuc_dt$end))

# Helper function to map methylation to tiles, calculate means, and add coverage labels
map_meth_to_tiles <- function(meth_df, method_name, tiles_gr, min_coverage) {
  
  # Filter by specified coverage
  filtered_df <- meth_df[meth_df$score >= min_coverage, ]
  
  # Convert to GRanges
  gr_meth <- GRanges(seqnames = filtered_df$chr,
                     ranges = IRanges(start = filtered_df$start, end = filtered_df$end),
                     beta = filtered_df$beta)
  
  # Find overlaps
  overlaps <- findOverlaps(gr_meth, tiles_gr)
  
  # Create mapping table
  meth_per_tile <- data.table(
    tile_id = subjectHits(overlaps),
    beta = filtered_df$beta[queryHits(overlaps)]
  )
  
  # Calculate mean per tile and add labels
  tile_means <- meth_per_tile[, .(mean_meth = mean(beta, na.rm = TRUE)), by = tile_id]
  tile_means[, Method := method_name]
  tile_means[, Coverage := paste0(min_coverage, "X")]
  
  return(tile_means)
}

# 2. Loop through paired samples
for (i in seq_along(file_paths_ONT)) {
  
  sample_id <- str_extract(basename(file_paths_ONT[i]), "\\d+")
  message(sprintf("Processing Sample: %s...", sample_id))
  
  # Load matched data into memory
  ont_data <- readRDS(file_paths_ONT[i])
  em_data  <- readRDS(file_paths_EM[i])
  
  message("  Mapping ONT data (10X and 1X)...")
  ont_10x <- map_meth_to_tiles(ont_data, "Nanopore", gr_tiles, 10)
  em_10x  <- map_meth_to_tiles(em_data, "EM-Seq", gr_tiles, 10)
  
  # Combine all the aggregated data
  combined_means <- rbind(ont_10x, em_10x)
  
  # Set factors to control facet layout order
  combined_means[, Coverage := factor(Coverage, levels = c("10X"))] 
  combined_means[, Method := factor(Method, levels = c("EM-Seq", "Nanopore"))]
  
  # Merge the calculated means back with the GC% data
  plot_data <- merge(nuc_dt, combined_means, by = "tile_id", all.y = TRUE, allow.cartesian = TRUE) 
  plot_data <- plot_data[!is.na(mean_meth)]
  
  # ==========================================
  # NEW: Create a dataframe for the A & B tags
  # ==========================================
  tag_df <- data.frame(
    Coverage = factor("10X", levels = c("10X")), # Matches your single coverage level
    Method = factor(c("EM-Seq", "Nanopore"), levels = c("EM-Seq", "Nanopore")),
    tag = c("A", "B")
  )
  
  # ==========================================
  # 3. Plotting: Methylation vs. GC Content
  # ==========================================
  message("  Generating plot...")
  
  m_gc <- ggplot(plot_data, aes(x = pct_gc, y = mean_meth)) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(trans = "log10", name = "Tile Count\n(log10)") +
    scale_y_continuous(limits = c(-0.01, 1.01)) + 
    scale_x_continuous(limits = c(-0.01, 1.01)) + 
    facet_grid(~ Method) +
    # Add the tags here!
    geom_text(
      data = tag_df, 
      aes(x = -Inf, y = Inf, label = tag), 
      hjust = -0.5,   # Push slightly right from the left edge
      vjust = 1.5,    # Push slightly down from the top edge
      size = 6,       # Adjust size as needed
      fontface = "bold", 
      inherit.aes = FALSE # Crucial so it doesn't look for pct_gc/mean_meth
    ) +
    labs(
      title = NULL,
      x = "GC Content (%)",
      y = "Average Methylation (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
  
  # Save the plot
  save_path <- file.path(out_dir, paste0("methyl_GC_EM_vs_ONT_", sample_id, ".png"))
  ggsave(save_path, plot = m_gc, width = 10, height = 6, dpi = 300)
  
  # Clean up memory
  rm(ont_data, em_data, ont_10x, em_10x, combined_means, plot_data, m_gc, tag_df)
  gc()
}

message("All samples processed successfully.")