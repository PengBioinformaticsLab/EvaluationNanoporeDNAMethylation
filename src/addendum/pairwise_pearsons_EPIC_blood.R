library(tidyverse)

out_dir <- "addendum/res"
if (!dir.exists(out_dir)) dir.create(out_dir)

# --- Blood: EPIC-only ---
# File Loading

sample_match <- data.frame(
  nanopore = c("s4059246541", "s4059246577", "s4059249475", "s4059249500"),
  epic = c("RID_4501", "RID_4077", "RID_4466", "RID_4390")
)



### Load Data

Epic_Manifest <- read.table(
  gzfile("../nanopore_epic/data/EPIC.hg38.manifest.tsv.gz"),
  header = TRUE, sep = ""
)



Beta_Epic <- read.csv(
  "../nanopore_epic/data/4_individuals_betas 2.csv"
)

cn <- colnames(Beta_Epic)
cn[1] <- "CpG"
colnames(Beta_Epic) <- cn

idx <- stringr::str_detect(Beta_Epic$CpG, pattern = "cg")
Beta_Epic <- Beta_Epic[idx, ]
Beta_Epic <- Beta_Epic[Beta_Epic$CpG %in% Epic_Manifest$Probe_ID, ]
Beta_Epic <- Beta_Epic[, c("CpG", sample_match$epic)]

idx <- match(Beta_Epic$CpG, Epic_Manifest$Probe_ID)
Epic_Manifest <- Epic_Manifest[idx, ]

# 1. Identify the sample columns from your matched metadata
sample_cols <- sample_match$epic

# 2. Get pairwise combinations of those column names
epic_pairs <- combn(sample_cols, 2, simplify = FALSE)

# 3. Execute comparisons and build the final table directly from columns
epic_tab <- map_dfr(epic_pairs, function(p) {
  samp1 <- p[1]
  samp2 <- p[2]

  # Extract the columns directly as vectors (no joining needed!)
  vec1 <- Beta_Epic[[samp1]]
  vec2 <- Beta_Epic[[samp2]]

  # Since there are no NAs and columns are aligned, the count is simply the row count
  n_cpgs <- nrow(Beta_Epic)

  # Compute Correlation and CIs
  if (n_cpgs > 2) {
    c_test <- suppressWarnings(cor.test(vec1, vec2, method = "pearson"))
    r_est <- round(c_test$estimate, 3)
    ci_low <- round(c_test$conf.int[1], 3)
    ci_high <- round(c_test$conf.int[2], 3)
  } else {
    r_est <- NA_real_
    ci_low <- NA_real_
    ci_high <- NA_real_
  }

  # Return tibble row
  tibble(
    comparison = paste(samp1, "vs.", samp2),
    cpgs_detected = n_cpgs, # Replaces X, Y, and Overlap since they are identical here
    pearson_cor = r_est,
    pearson_ci_lower = ci_low,
    pearson_ci_upper = ci_high
  )
})

# 4. Add platform metadata and save
epic_tab <- epic_tab %>% mutate(platform = "EPIC")

write_csv(epic_tab, file.path(out_dir, "epic_blood_paired_pearsons.csv"))

rm(epic_tab)
gc()
