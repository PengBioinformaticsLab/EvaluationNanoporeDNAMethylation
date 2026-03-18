library(tidyverse)
library(hexbin)
library(patchwork)

out_dir <- "addendum/res"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ==========================================
# 1. Data Loading and Prep
# ==========================================

# --- Brain Data ---
nano_df_brain <- readRDS("data/nanopore/nano_51.rds") |>
  mutate(pos = paste(chr, start + 1, sep = "_"))

em_df <- readRDS("data/em-seq/em_51.rds") |>
  mutate(pos = paste(chr, start, sep = "_"))

# --- Blood Data ---
nano_df_blood <- readRDS("../nanopore_epic/data/sample_4059246577_analysis/s4059246577_combined.rds") |>
  mutate(pos = paste(chr, start, sep = "_"))

sample_match <- data.frame(
  nanopore = c("s4059246541", "s4059246577", "s4059249475", "s4059249500"),
  epic = c("RID_4501", "RID_4077", "RID_4466", "RID_4390")
)

Epic_Manifest <- read.table(
  gzfile("../nanopore_epic/data/EPIC.hg38.manifest.tsv.gz"),
  header = TRUE, sep = ""
)

Beta_Epic <- read.csv("../nanopore_epic/data/4_individuals_betas 2.csv")

cn <- colnames(Beta_Epic)
cn[1] <- "CpG"
colnames(Beta_Epic) <- cn

idx <- stringr::str_detect(Beta_Epic$CpG, pattern = "cg")
Beta_Epic <- Beta_Epic[idx, ]
Beta_Epic <- Beta_Epic[Beta_Epic$CpG %in% Epic_Manifest$Probe_ID, ]
Beta_Epic <- Beta_Epic[, c("CpG", sample_match$epic)]

idx <- match(Beta_Epic$CpG, Epic_Manifest$Probe_ID)
Epic_Manifest <- Epic_Manifest[idx, ]

Beta_Epic <- Beta_Epic[, c("CpG", "RID_4077")]
Epic_Manifest$pos <- paste(Epic_Manifest$CpG_chrm, Epic_Manifest$CpG_beg, sep = "_")

epic_df <- Beta_Epic %>%
  left_join(Epic_Manifest %>% dplyr::select(Probe_ID, pos),
    by = c("CpG" = "Probe_ID")
  ) %>%
  transmute(pos, beta = RID_4077)


# ==========================================
# 2. Execution, Combining, and Plotting
# ==========================================
cutoff <- 10

# --- Process Brain ---
nano_filt_brain <- nano_df_brain %>%
  filter(score >= cutoff) %>%
  mutate(beta = (NMod + NOtherMod) / score)

em_filt <- em_df %>% filter(score >= cutoff)

shared_brain <- inner_join(
  nano_filt_brain %>% select(pos, beta_nano = beta),
  em_filt %>% select(pos, beta_ref = beta),
  by = "pos"
) %>%
  mutate(mean = (beta_nano + beta_ref) / 2, diff = beta_nano - beta_ref)

stats_brain <- shared_brain %>%
  summarize(
    mean_diff = mean(diff, na.rm = TRUE),
    sd_diff   = sd(diff, na.rm = TRUE),
    loa_upper = mean_diff + 1.96 * sd_diff,
    loa_lower = mean_diff - 1.96 * sd_diff
  )

# --- Plot Brain ---
plot_brain <- ggplot(shared_brain, aes(x = mean, y = diff)) +
  stat_binhex(aes(fill = after_stat(ndensity)), bins = 300) +
  scale_fill_viridis_c(
    option = "mako", 
    direction = -1, 
    trans = "log10", 
    breaks = function(x) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)), 
    labels = c("Min", "Max") 
  ) +
  geom_hline(yintercept = stats_brain$mean_diff, color = "#FF5722", linewidth = 1, alpha = 0.3) +
  geom_hline(yintercept = stats_brain$loa_upper, color = "#FF5722", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = stats_brain$loa_lower, color = "#FF5722", linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(limits = c(-1.01, 1.01)) + 
  theme_minimal() +
  labs(
    title = "",
    x = "(Nanopore + EM-Seq) / 2",
    y = "Nanopore - EM-Seq",
    fill = "Density"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# --- Process Blood ---
nano_filt_blood <- nano_df_blood %>%
  filter(score >= cutoff) %>%
  mutate(beta = (NMod + NOtherMod) / score)

shared_blood <- inner_join(
  nano_filt_blood %>% select(pos, beta_nano = beta),
  epic_df %>% select(pos, beta_ref = beta),
  by = "pos"
) %>%
  mutate(mean = (beta_nano + beta_ref) / 2, diff = beta_nano - beta_ref)


stats_blood <- shared_blood %>%
  summarize(
    mean_diff = mean(diff, na.rm = TRUE),
    sd_diff   = sd(diff, na.rm = TRUE),
    loa_upper = mean_diff + 1.96 * sd_diff,
    loa_lower = mean_diff - 1.96 * sd_diff
  )

# --- Plot Blood ---
plot_blood <- ggplot(shared_blood, aes(x = mean, y = diff)) +
  stat_binhex(aes(fill = after_stat(ndensity)), bins = 300) +
  scale_fill_viridis_c(
    option = "rocket", 
    direction = -1, 
    trans = "log10", 
    breaks = function(x) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)), 
    labels = c("Min", "Max") 
  ) +
  geom_hline(yintercept = stats_blood$mean_diff, color = "#22A8FF", linewidth = 1, alpha = 0.3) +
  geom_hline(yintercept = stats_blood$loa_upper, color = "#22A8FF", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = stats_blood$loa_lower, color = "#22A8FF", linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(limits = c(-1.01, 1.01)) + 
  theme_minimal() +
  labs(
    title = "",
    x = "(Nanopore + EPIC) / 2",
    y = "Nanopore - EPIC",
    fill = "Density"
  ) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# --- Combine with Patchwork ---
combined_plot <- ( plot_brain + plot_blood) +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = "bottom")


# --- Save ---
ggsave(
  file.path(out_dir, paste0("Bland-Altman_Combined_cov", cutoff, "X.png")),
  plot = combined_plot,
  width = 12,
  height = 6
)
# 