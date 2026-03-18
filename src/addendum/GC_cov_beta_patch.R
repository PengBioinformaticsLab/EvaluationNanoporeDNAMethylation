# Load libraries
library(data.table)
library(GenomicRanges)
library(ggplot2)
library(stringr)
library(patchwork) 

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
nuc_dt[, tile_id := .I]

gr_tiles <- GRanges(seqnames = nuc_dt$chr,
                    ranges = IRanges(start = nuc_dt$start, end = nuc_dt$end))

# --- Helper Functions ---

# Helper 1: Coverage
map_cov_to_tiles <- function(df, platform_name, tiles_gr) {
  mean_score <- mean(df$score, na.rm = TRUE)
  df$norm_coverage <- df$score / mean_score
  
  gr_cov <- GRanges(seqnames = df$chr,
                    ranges = IRanges(start = df$start, end = df$start),
                    norm_coverage = df$norm_coverage)
  
  overlaps <- findOverlaps(gr_cov, tiles_gr)
  cov_per_tile <- data.table(
    tile_id = subjectHits(overlaps),
    norm_coverage = gr_cov$norm_coverage[queryHits(overlaps)]
  )
  
  tile_means <- cov_per_tile[, .(mean_norm_cov = mean(norm_coverage, na.rm = TRUE)), by = tile_id]
  tile_means[, Platform := platform_name]
  return(tile_means)
}

# Helper 2: Methylation 
map_meth_to_tiles <- function(meth_df, tiles_gr, min_coverage) {
  filtered_df <- meth_df[meth_df$score >= min_coverage, ]
  
  gr_meth <- GRanges(seqnames = filtered_df$chr,
                     ranges = IRanges(start = filtered_df$start, end = filtered_df$end),
                     beta = filtered_df$beta)
  
  overlaps <- findOverlaps(gr_meth, tiles_gr)
  meth_per_tile <- data.table(
    tile_id = subjectHits(overlaps),
    beta = filtered_df$beta[queryHits(overlaps)]
  )
  
  tile_means <- meth_per_tile[, .(mean_meth = mean(beta, na.rm = TRUE)), by = tile_id]
  tile_means[, Coverage := paste0(min_coverage, "X")]
  return(tile_means)
}

# 2. Loop through paired samples
for (i in seq_along(file_paths_ONT)) {
  
  sample_id <- str_extract(basename(file_paths_ONT[i]), "\\d+")
  message(sprintf("Processing Sample: %s...", sample_id))
  
  ont_data <- readRDS(file_paths_ONT[i])
  em_data  <- readRDS(file_paths_EM[i])
  
  # Map Coverage 
  message("  Mapping & Normalizing Coverage...")
  ont_cov <- map_cov_to_tiles(ont_data, "Nanopore", gr_tiles)
  em_cov  <- map_cov_to_tiles(em_data, "EM-Seq", gr_tiles)
  
  # Map Methylation (1X and 10X)
  message("  Mapping Methylation (1X and 10X)...")
  ont_meth_1x  <- map_meth_to_tiles(ont_data, gr_tiles, 1)
  em_meth_1x   <- map_meth_to_tiles(em_data, gr_tiles, 1)
  
  ont_meth_10x <- map_meth_to_tiles(ont_data, gr_tiles, 10)
  em_meth_10x  <- map_meth_to_tiles(em_data, gr_tiles, 10)
  
  # Merge Data
  message("  Merging data tables...")
  
  # 1X Data Merge
  ont_comb_1x <- merge(ont_cov, ont_meth_1x, by = "tile_id", all.x = TRUE)
  em_comb_1x  <- merge(em_cov, em_meth_1x, by = "tile_id", all.x = TRUE)
  data_1x <- rbind(ont_comb_1x, em_comb_1x)
  data_1x[, Platform := factor(Platform, levels = c("EM-Seq", "Nanopore"))]
  plot_data_1x <- merge(nuc_dt, data_1x, by = "tile_id", all.y = TRUE)
  
  # 10X Data Merge
  ont_comb_10x <- merge(ont_cov, ont_meth_10x, by = "tile_id", all.x = TRUE)
  em_comb_10x  <- merge(em_cov, em_meth_10x, by = "tile_id", all.x = TRUE)
  data_10x <- rbind(ont_comb_10x, em_comb_10x)
  data_10x[, Platform := factor(Platform, levels = c("EM-Seq", "Nanopore"))]
  plot_data_10x <- merge(nuc_dt, data_10x, by = "tile_id", all.y = TRUE)
  
  # ==========================================
  # 3. Plotting using Patchwork
  # ==========================================
  message("  Generating combined plots...")
  
  # Calculate universal 99.9th percentile to cap Y-axis (shared across all)
  plot_data_cov <- plot_data_1x[!is.na(mean_norm_cov)] # Cov is identical for 1X/10X
  y_cap <- quantile(plot_data_cov$mean_norm_cov, 0.999, na.rm = TRUE)
  
  # Plot 1: Top Row (Tile Count log10 - Generated once per sample)
  p_top <- ggplot(plot_data_cov, aes(x = pct_gc, y = mean_norm_cov)) +
    geom_hex(bins = 100) +
    scale_fill_viridis_c(trans = "log10", name = "Frequency\n(log10)") +
    scale_y_continuous(limits = c(0, y_cap)) + 
    scale_x_continuous(limits = c(-0.01, 1.01)) + 
    facet_wrap(~ Platform) +
    labs(x = NULL, y = "Mean Norm. Coverage") +
    theme_minimal() +
    theme(
      axis.title.y = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
  
  
  
  # Helper function to generate the bottom methylation plot
  generate_meth_plot <- function(df) {
    ggplot(df[!is.na(mean_norm_cov) & !is.na(mean_meth)], aes(x = pct_gc, y = mean_norm_cov, z = mean_meth)) +
      stat_summary_hex(fun = mean, bins = 100) +
      scale_fill_viridis_c(option = "magma", name = "Mean\nMethylation") +
      scale_y_continuous(limits = c(0, y_cap)) + 
      scale_x_continuous(limits = c(-0.01, 1.01)) + 
      facet_wrap(~ Platform) +
      labs(x = "GC Content (%)", y = "Mean Norm. Coverage") +
      theme_minimal() +
      theme(
        axis.title = element_text(face = "bold"),
        strip.text = element_blank(), 
        legend.position = "right"
      )
  }
  
  # Generate bottom plots
  p_bottom_1x  <- generate_meth_plot(plot_data_1x)
  
  # Combine and Save 1X Plot
  p_combined_1x <- p_top / p_bottom_1x + 
    plot_annotation(
      tag_levels = 'A',  
      title = NULL,
      subtitle = NULL,
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  ggsave(file.path(out_dir, paste0("Combined_GC_Cov_Meth_1X_", sample_id, ".png")), 
         plot = p_combined_1x, width = 11, height = 9, dpi = 300)
  
  
  # Clean up memory
  rm(ont_data, em_data, ont_cov, em_cov, 
     ont_meth_1x, em_meth_1x, ont_meth_10x, em_meth_10x,
     ont_comb_1x, em_comb_1x, ont_comb_10x, em_comb_10x, 
     data_1x, data_10x, plot_data_1x, plot_data_10x, plot_data_cov, 
     p_top, p_bottom_1x, p_bottom_10x, p_combined_1x, p_combined_10x)
  gc()
}

message("All samples processed successfully.")