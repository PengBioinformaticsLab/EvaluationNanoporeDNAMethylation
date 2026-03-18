library(tidyverse)
library(openxlsx)

out_dir <- "addendum/res"
if (!dir.exists(out_dir)) dir.create(out_dir)
     

EM_brain <- read.csv("addendum/res/EM_brain_paired_pearsons.csv")
EM_brain$tissue <- "Brain"

Nano_brain <- read.csv("addendum/res/nanopore_brain_paired_pearsons.csv")
Nano_brain$tissue <- "Brain"

Nano_blood <- read.csv("addendum/res/nanopore_blood_paired_pearsons.csv")
Nano_blood$tissue <- "Blood"

sum_df <- rbind(EM_brain, Nano_brain, Nano_blood) |> 
  mutate(cutoff_coverage = as.character(cutoff_coverage))

EPIC_blood <- read.csv("addendum/res/epic_blood_paired_pearsons.csv")
EPIC_blood$tissue <- "Blood"


EPIC_blood <- EPIC_blood |> 
  mutate(X_cpgs_detected = cpgs_detected,
         Y_cpgs_detected = cpgs_detected, 
         overlap_cpgs_detected = cpgs_detected,
         cutoff_coverage = "Array (EPIC)") |> 
  dplyr::select(-cpgs_detected)



# make cols missing in EPIC NA and join all of the above together 
paired_df <- bind_rows(sum_df, EPIC_blood)


write.csv(paired_df, "addendum/res/intersample_pearsons.csv")


paired_df <- paired_df |> 
  transmute(
    pearson_cor, 
    cutoff_coverage = as.character(cutoff_coverage), 
    tissue, 
    platform, 
    comparison, 
    study = "paired"
  )




# Original platform vs. platform comparison

og_brain <- read.csv("output/tables/cpg_stats.csv") |> 
  janitor::clean_names() |> 
  transmute(cutoff_coverage = coverage_threshold,
            pearson_cor = correlation_pearson,
            tissue = "Brain",
            platform = "Nanopore * EM-Seq",
            comparison = sample,
            study = "original") |> 
  mutate(cutoff_coverage = as.character(parse_number(cutoff_coverage)))


og_blood <- read.csv("../nanopore_epic/rlt/tables/cor_beta.csv") |> 
  transmute(cutoff_coverage = "Array (EPIC)",
            pearson_cor = corr,
            tissue = "Blood",
            platform = "Nanopore * EPIC",
            comparison = X,
            study = "original")

combined_df <- bind_rows(paired_df, og_blood, og_brain) |>
  mutate(
    cutoff_coverage = as.character(cutoff_coverage),
    cutoff_coverage = if_else(is.na(cutoff_coverage) | cutoff_coverage == "NA", 
                              "Array (EPIC)", 
                              cutoff_coverage),
    # Now it's safe to factor
    cutoff_coverage = factor(cutoff_coverage, levels = c("1", "5", "10", "15", "20","Array (EPIC)"))
  ) 



# Using unsummarized raw data
paired_cor_plot <- combined_df |>
  # 2. Setup aesthetics mapping
  ggplot(aes(
    x = cutoff_coverage, 
    y = pearson_cor
  )) +
  # 4. Add the jittered points
  geom_jitter(
    aes(color = platform, shape = study),
    width = 0.2, # Horizontal spread
    height = 0,  # CRITICAL: Do not jitter the actual correlation values!
    size = 2.5,
    alpha = 0.7
  ) +
  # 5. Facet by tissue
  facet_wrap(~ tissue) +
  # 6. Styling and Labels
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) + 
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  ) +
  labs(
    title = "Correlation coefficients by coverage, platform, and tissue",
    x = "Cutoff Coverage",
    y = "Pearson Correlation (r)",
    color = "",
    shape = "Analysis context"
  ) +
  scale_shape_manual(
    values = c("paired" = 16, "original" = 17), # 16 is a circle, 17 is a triangle
    labels = c("paired" = "Paired Replicates", "original" = "Cross-Platform Benchmark")
  )

ggsave(filename = file.path(out_dir, "paired_pearsons_plot.png"),
       plot = paired_cor_plot)


unique(combined_df$cutoff_coverage)

# GP table
# For {tissue}
for (t in c("Brain")){
    
  c <- 10
    
    df <- combined_df |> filter(cutoff_coverage == c & tissue == t)
    
    # Assuming your data frame is named 'df'
    
    # 1. Identify all unique sample names to determine the matrix dimensions
    # We use regex to find all occurrences of "Sample X" in the comparison column
    all_samples <- regmatches(df$comparison, gregexpr("Sample \\d+", df$comparison))
    samples <- sort(unique(unlist(all_samples)))
    
    # 2. Initialize an empty character matrix 
    cor_matrix <- matrix(NA_character_, nrow = length(samples), ncol = length(samples),
                         dimnames = list(samples, samples))
    
    # 3. Loop through each row of the dataframe and populate the matrix
    for (i in seq_len(nrow(df))) {
      platform <- trimws(df$platform[i])
      comp <- trimws(df$comparison[i])
      
      # Build the formatted string
      r_val <- df$pearson_cor[i]
      

        val <- sprintf("%.2f", r_val) 
      
      if (platform == "Nanopore * EM-Seq") {
        # DIAGONAL
        samp <- sub(" original", "", comp)
        cor_matrix[samp, samp] <- val
        
      } else {
        # OFF-DIAGONAL
        samps <- unlist(strsplit(comp, " vs. "))
        s1 <- samps[1]
        s2 <- samps[2]
        
        if (platform == "EM-Seq") {
          cor_matrix[s2, s1] <- val
        } else if (platform == "Nanopore") {
          cor_matrix[s1, s2] <- val
        }
      }
    }
    
    clean_names <- gsub(" ", "_", rownames(cor_matrix))
    
    # Do NOT convert to numeric here anymore. Just assign it directly.
    brain_cor_matrix <- cor_matrix
    
    rownames(brain_cor_matrix) <- paste0(clean_names)
    colnames(brain_cor_matrix) <- paste0(clean_names)
    print(brain_cor_matrix)
 }

#  Sample map

sample_match <- data.frame(
  epic = c("RID_4501", "RID_4077", "RID_4466", "RID_4390"),
  study = c("Sample_A", "Sample_B", "Sample_C", "Sample_D")
)


rid_mapping <- setNames(sample_match$study, sample_match$epic)

blood_df <- combined_df |>
  filter(cutoff_coverage %in% c(10, "Array (EPIC)") & tissue == "Blood")  |> 
  mutate(comparison = str_replace_all(comparison, rid_mapping))

# 1. Standardize the comparison strings to clean up prefixes and underscores
blood_df$clean_comp <- blood_df$comparison
blood_df$clean_comp <- gsub("Nanopore |EPIC |cor_beta_| original", "", blood_df$clean_comp)
blood_df$clean_comp <- gsub("Sample_", "Sample ", blood_df$clean_comp)

# 2. Identify all unique sample names
all_samples <- unlist(strsplit(blood_df$clean_comp, " vs. "))
samples <- sort(unique(all_samples))

# 3. Initialize an empty character matrix
cor_matrix <- matrix(NA_character_, nrow = length(samples), ncol = length(samples),
                     dimnames = list(samples, samples))

# 4. Loop through each row and populate the matrix
for (i in seq_len(nrow(blood_df))) {
  platform <- trimws(blood_df$platform[i])
  comp <- trimws(blood_df$clean_comp[i])
  
  # Build the formatted string
  r_val <- blood_df$pearson_cor[i]
  
  val <- sprintf("%.2f", r_val) 
  
  if (platform == "Nanopore * EPIC") {
    cor_matrix[comp, comp] <- val
  } else {
    samps <- unlist(strsplit(comp, " vs. "))
    s1 <- samps[1]
    s2 <- samps[2]
    
    if (s1 > s2) {
      temp <- s1; s1 <- s2; s2 <- temp
    }
    
    if (platform == "EPIC") {
      cor_matrix[s2, s1] <- val
    } else if (platform == "Nanopore") {
      cor_matrix[s1, s2] <- val
    }
  }
}

# Do NOT convert to numeric
blood_cor_matrix <- cor_matrix

clean_names <- gsub("Sample ", "Sample_", rownames(cor_matrix))
rownames(blood_cor_matrix) <- paste0(clean_names)
colnames(blood_cor_matrix) <- paste0(clean_names)
# View the updated matrix
print(blood_cor_matrix)


# Intersample_cor table ---------------------------------------------------

# --- 1. Prepare the DataFrames ---
# Convert Brain matrix
df_brain <- as.data.frame(brain_cor_matrix)
df_brain <- cbind("Reference Sample" = rownames(df_brain), df_brain)
rownames(df_brain) <- NULL
# Clean up column names (e.g., ONT_Sample_1 to ONT Sample 1)
colnames(df_brain) <- gsub("_", " ", colnames(df_brain))

# Convert Blood matrix
df_blood <- as.data.frame(blood_cor_matrix)
df_blood <- cbind("Reference Sample" = rownames(df_blood), df_blood)
rownames(df_blood) <- NULL
# Clean up column names
colnames(df_blood) <- gsub("_", " ", colnames(df_blood))

# --- 2. Create the Excel Workbook ---
wb <- createWorkbook()
addWorksheet(wb, "Table S1")

# Define styles for formatting
header_style <- createStyle(textDecoration = "bold", border = "bottom")
title_style <- createStyle(textDecoration = "bold", fontSize = 12)

# NEW: Define specific styles for the matrix regions
bold_style <- createStyle(textDecoration = "bold", numFmt = "0.00")      # Diagonal (Cross-platform)
italic_style <- createStyle(textDecoration = "italic", numFmt = "0.00")  # Lower left (Orthogonal)
standard_num_style <- createStyle(numFmt = "0.00")                       # Upper right (Nanopore)

# Define a style for the long caption (Wrap text, align top left)
caption_style <- createStyle(wrapText = TRUE, halign = "left", valign = "top")

# --- 2.5 Add the Caption ---
# (Updated caption to include the formatting legend)
caption_text <- "Supplementary Table S1: Cross-platform and within-platform Pearson correlations of DNA methylation. Matrices display Pearson correlation coefficients (r) for matched samples in Brain (Panel A) and Blood (Panel B) at a cutoff coverage of 10x. For both matrices: the main diagonal (bold) represents the cross-platform correlation for the exact same sample (e.g., EM-seq vs. Nanopore). The upper right triangle (standard) displays the within-platform correlations for Nanopore between different samples. The lower left triangle (italics) displays the within-platform correlations for the orthogonal technology (EM-seq in Panel A; EPIC Array in Panel B)."

# Write the caption to the top-left cell
writeData(wb, "Table S1", caption_text, startRow = 1, startCol = 1)

# Merge cells across Row 1 to match the width of your dataframe
mergeCells(wb, "Table S1", cols = 1:ncol(df_brain), rows = 1)

# Apply the wrapping style and manually set a taller row height so it fits
addStyle(wb, "Table S1", caption_style, rows = 1, cols = 1)
setRowHeights(wb, "Table S1", rows = 1, heights = 110) 


# --- 3. Write Brain Data (Panel A) ---
start_a_row <- 3

# Add Panel A Title
writeData(wb, "Table S1", "A. Brain (EM-seq vs. Nanopore)", startRow = start_a_row, startCol = 1)
addStyle(wb, "Table S1", title_style, rows = start_a_row, cols = 1)

# Add Brain DataFrame
writeData(wb, "Table S1", df_brain, startRow = start_a_row + 1, startCol = 1, headerStyle = header_style)

# Apply Bold/Italic formatting to Panel A
n_samples_brain <- nrow(df_brain)
for (i in 1:n_samples_brain) {
  for (j in 1:n_samples_brain) {
    excel_r <- start_a_row + 1 + i  # +1 offsets the header row
    excel_c <- 1 + j                # +1 offsets the 'Reference Sample' column
    
    if (i == j) {
      addStyle(wb, "Table S1", bold_style, rows = excel_r, cols = excel_c)
    } else if (i > j) {
      addStyle(wb, "Table S1", italic_style, rows = excel_r, cols = excel_c)
    } else {
      addStyle(wb, "Table S1", standard_num_style, rows = excel_r, cols = excel_c)
    }
  }
}

# --- 4. Write Blood Data (Panel B) ---
# Calculate where to start Panel B dynamically
start_b_row <- start_a_row + 1 + nrow(df_brain) + 2

# Add Panel B Title
writeData(wb, "Table S1", "B. Blood (EPIC vs. Nanopore)", startRow = start_b_row, startCol = 1)
addStyle(wb, "Table S1", title_style, rows = start_b_row, cols = 1)

# Add Blood DataFrame
writeData(wb, "Table S1", df_blood, startRow = start_b_row + 1, startCol = 1, headerStyle = header_style)

# Apply Bold/Italic formatting to Panel B
n_samples_blood <- nrow(df_blood)
for (i in 1:n_samples_blood) {
  for (j in 1:n_samples_blood) {
    excel_r <- start_b_row + 1 + i
    excel_c <- 1 + j
    
    if (i == j) {
      addStyle(wb, "Table S1", bold_style, rows = excel_r, cols = excel_c)
    } else if (i > j) {
      addStyle(wb, "Table S1", italic_style, rows = excel_r, cols = excel_c)
    } else {
      addStyle(wb, "Table S1", standard_num_style, rows = excel_r, cols = excel_c)
    }
  }
}

# --- 5. Polish and Save ---
# Auto-adjust column widths so the sample names aren't cut off
setColWidths(wb, "Table S1", cols = 1:ncol(df_brain), widths = "auto")

# Bold the first column (Reference Samples) for both matrices to make it look clean
addStyle(wb, "Table S1", createStyle(textDecoration = "bold"), 
         rows = (start_a_row + 2):(start_a_row + 1 + n_samples_brain), cols = 1)
addStyle(wb, "Table S1", createStyle(textDecoration = "bold"), 
         rows = (start_b_row + 2):(start_b_row + 1 + n_samples_blood), cols = 1)

# Save the file
saveWorkbook(wb, "addendum/res/Supplementary_Table_S1.xlsx", overwrite = TRUE)