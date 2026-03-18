library(tidyverse)

out_dir <- "addendum/res"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Helper Function
# Compute correlation and CI between sample X and sample Y across coverage cutoffs
count_cpgs <- function(samp_X, samp_Y) {
  
  cutoff_coverage <- c(1, 5, 10, 15, 20)
  
  # map_dfr loops through cutoffs and returns a tidy data frame automatically
  res_df <- map_dfr(cutoff_coverage, function(cutoff) {
    
    # 1. Filter by cutoff
    X_filt <- samp_X %>% filter(score >= cutoff)
    Y_filt <- samp_Y %>% filter(score >= cutoff)
    
    # 2. Inner join ensures perfect alignment of 'pos' before correlation
    shared <- inner_join(
      X_filt %>% select(pos, beta_X = beta),
      Y_filt %>% select(pos, beta_Y = beta),
      by = "pos"
    )
    
    # 3. Compute Correlation and CIs (Needs >2 observations for CI)
    if (nrow(shared) > 2) {
      # suppressWarnings handles edge cases where beta values have zero variance
      c_test <- suppressWarnings(cor.test(shared$beta_X, shared$beta_Y, method = "pearson"))
      r_est   <- round(c_test$estimate, 3)
      ci_low  <- round(c_test$conf.int[1], 3)
      ci_high <- round(c_test$conf.int[2], 3)
    } else {
      r_est   <- NA_real_
      ci_low  <- NA_real_
      ci_high <- NA_real_
    }
    
    # 4. Return a one-row tibble for this specific cutoff
    tibble(
      cutoff_coverage = cutoff,
      X_cpgs_detected = nrow(X_filt),
      Y_cpgs_detected = nrow(Y_filt),
      overlap_cpgs_detected = nrow(shared),
      pearson_cor = r_est,
      pearson_ci_lower = ci_low,
      pearson_ci_upper = ci_high
    )
  })
  
  return(res_df)
}

# --- Brain: Nanopore-only ---

# File Loading
file_paths <- list.files(path = "data/nanopore", pattern = "nano_.*\\.rds$", full.names = TRUE)
file_order <- order(as.numeric(str_extract(basename(file_paths), "\\d+")))
file_paths <- file_paths[file_order]

# Read all files into a list and name them automatically (e.g., Sample_1, Sample_2)
sample_names <- paste("Sample", seq_along(file_paths))
samples <- map(file_paths, readRDS) %>% set_names(sample_names)

# Pre-process Data (Compute 'pos' once per sample to save time)
samples <- map(samples, ~ mutate(.x, pos = paste(chr, start, sep = "_")))

# Get pairwise combinations
pairs <- combn(names(samples), 2, simplify = FALSE)

# Execute comparisons and build the final table
# map_dfr will run the comparisons and immediately bind them into one tall data frame
nano_tab <- map_dfr(pairs, function(p) {
  p1 <- p[1]
  p2 <- p[2]
  
  # Run function and add comparison column
  res <- count_cpgs(samples[[p1]], samples[[p2]])
  res <- res %>% mutate(comparison = paste(p1, "vs.", p2))
  
  return(res)
})

# Add platform metadata and save
nano_tab <- nano_tab %>% mutate(platform = "Nanopore")

write_csv(nano_tab, file.path(out_dir, "nanopore_brain_paired_pearsons.csv"))

rm(nano_tab)
gc()


# --- Brain: EM-only ---

# File Loading
file_paths <- list.files(path = "data/em-seq", pattern = "em_.*\\.rds$", full.names = TRUE)
file_order <- order(as.numeric(str_extract(basename(file_paths), "\\d+")))
file_paths <- file_paths[file_order]

# Read all files into a list and name them automatically (e.g., Sample_1, Sample_2)
sample_names <- paste("Sample", seq_along(file_paths))
samples <- map(file_paths, readRDS) %>% set_names(sample_names)

# Pre-process Data (Compute 'pos' once per sample to save time)
samples <- map(samples, ~ mutate(.x, pos = paste(chr, start, sep = "_")))

# Get pairwise combinations
pairs <- combn(names(samples), 2, simplify = FALSE)

# Execute comparisons and build the final table
# map_dfr will run the comparisons and immediately bind them into one tall data frame
em_tab <- map_dfr(pairs, function(p) {
  p1 <- p[1]
  p2 <- p[2]
  
  # Run function and add comparison column
  res <- count_cpgs(samples[[p1]], samples[[p2]])
  res <- res %>% mutate(comparison = paste(p1, "vs.", p2))
  
  return(res)
})

# Add platform metadata and save
em_tab <- em_tab %>% mutate(platform = "EM-Seq")

write_csv(em_tab, file.path(out_dir, "EM_brain_paired_pearsons.csv"))
