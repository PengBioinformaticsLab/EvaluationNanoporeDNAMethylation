library(openxlsx)
library(dplyr)
library(tidyr)

# Input / Output Setup
out_dir <- "addendum/res"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ==========================================
# 1. Load Data
# ==========================================
# EM-Seq * ONT
supp_table_EM <- read.csv(file.path(out_dir, "positive_agreement_summary_EM.csv"))
supp_table_EM$Sample <- paste0("Sample ", as.numeric(factor(supp_table_EM$Sample_ID)))
supp_table_EM$Sample_ID <- NULL
conf_mat_EM <- readRDS(file.path(out_dir, "confusion_matrices_EM.rds"))

# EPIC * ONT
epic_dir <- "../nanopore_epic/rlt/tables/"
supp_table_EPIC <- read.csv(file.path(epic_dir, "binned_methyl_agreement_EPIC.csv"))
conf_mat_EPIC <- readRDS(file.path(epic_dir, "confusion_matrices_EPIC.rds"))

# Rename the required EM samples in the list
names(conf_mat_EM)[names(conf_mat_EM) == "51"] <- "Sample 1"
names(conf_mat_EM)[names(conf_mat_EM) == "69"] <- "Sample 2"
names(conf_mat_EM)[names(conf_mat_EM) == "168"] <- "Sample 3" 
names(conf_mat_EM)[names(conf_mat_EM) == "266"] <- "Sample 4"

# ==========================================
# 2. Format Confusion Matrices
# ==========================================
# Helper function to maintain the 3x3 layout
format_3x3_matrices <- function(mat_list, dataset_name) {
  do.call(rbind, lapply(names(mat_list), function(samp) {
    mat <- mat_list[[samp]]
    
    data.frame(
      Dataset = dataset_name,
      Sample = samp,
      Assay_Level = c("Low", "Medium", "High"), 
      ONT_Low = mat[, "Low"],                   
      ONT_Medium = mat[, "Medium"],             
      ONT_High = mat[, "High"],                 
      row.names = NULL
    )
  }))
}

# Process both datasets and stack into one final table
df_EM <- format_3x3_matrices(conf_mat_EM, "ONT vs. EM")
df_EPIC <- format_3x3_matrices(conf_mat_EPIC, "ONT vs. EPIC")
combined_conf_mat <- rbind(df_EM, df_EPIC) 



# ==========================================
# 3. Calculate Precision (Row Percentages)
# ==========================================
#Precision = (Matching ONT Calls) / (Total Assay Calls for that Bin) # PPA or "Class Agreement"
precision_df <- combined_conf_mat %>%
  mutate(
    Total_Assay = ONT_Low + ONT_Medium + ONT_High,
    Precision = case_when(
      Assay_Level == "Low" ~ ONT_Low / Total_Assay,
      Assay_Level == "Medium" ~ ONT_Medium / Total_Assay,
      Assay_Level == "High" ~ ONT_High / Total_Assay
    )
  ) %>%
  select(Dataset, Sample, Assay_Level, Precision) %>%
  pivot_wider(names_from = Assay_Level, values_from = Precision, names_prefix = "Precision_")



# ==========================================
# 5. Build Combined Unified Table
# ==========================================
supp_table_EM$Dataset <- "ONT vs. EM"
supp_table_EPIC$Dataset <- "ONT vs. EPIC"
combined_table_base <- rbind(supp_table_EM, supp_table_EPIC) |> 
  dplyr::select(-starts_with("PA_"))

# Join metrics and pivot them into long format (Bin: Low, Medium, High)
metrics_df <- combined_table_base %>%
  inner_join(precision_df, by = c("Dataset", "Sample")) %>%
  pivot_longer(
    cols = contains("_Low") | contains("_Medium") | contains("_High"),
    names_to = c(".value", "Bin"),
    names_sep = "_"
  )

# Merge Confusion Matrix counts with Metrics
final_merged <- combined_conf_mat %>%
  inner_join(metrics_df, by = c("Dataset", "Sample", "Assay_Level" = "Bin"))

# Construct the final display format row-by-row
final_display_rows <- list()

for(assay_name in unique(final_merged$Dataset)) {
  # Strip the prefix to get just the assay name (e.g., "EM" or "EPIC")
  short_assay <- gsub("ONT vs. ", "", assay_name)
  
  # 1. Add Assay Header Row (Corrected to 5 columns since Class Agreement is merged)
  final_display_rows[[length(final_display_rows) + 1]] <- data.frame(
    `Assay & Sample` = short_assay,
    `Assay Level` = NA,
    `ONT Low` = NA, 
    `ONT Inter` = NA, 
    `ONT High` = NA,
    check.names = FALSE
  )
  
  assay_data <- final_merged %>% filter(Dataset == assay_name)
  
  # 2. Add Sample Rows
  for(samp in unique(assay_data$Sample)) {
    samp_data <- assay_data %>% filter(Sample == samp)
    
    # Ensure correct order of the bins
    samp_data <- samp_data %>%
      mutate(Assay_Level = factor(Assay_Level, levels = c("Low", "Medium", "High"))) %>%
      arrange(Assay_Level)
    
    # Format the combined Sample + Agreement string for the first row of each sample block
    overall_acc <- round(samp_data$Overall_Accuracy[1], 3)
    samp_label <- paste0(samp, "\n(Agreement : ", overall_acc, ")")
    
    for(i in 1:nrow(samp_data)) {
      lvl <- as.character(samp_data$Assay_Level[i])
      display_lvl <- ifelse(lvl == "Medium", "Intermediate", lvl)
      
      # Combine Assay Level and Class Agreement (Precision) into one string
      class_agree_val <- sprintf("%.3f", samp_data$Precision[i])
      combined_assay_lvl <- paste0(display_lvl, "\n(CA: ", class_agree_val, ")")
      
      row_df <- data.frame(
        `Assay & Sample` = ifelse(i == 1, samp_label, NA), # Only print name on the first row
        `Assay Level` = combined_assay_lvl,
        `ONT Low` = samp_data$ONT_Low[i],
        `ONT Inter` = samp_data$ONT_Medium[i],
        `ONT High` = samp_data$ONT_High[i],
        check.names = FALSE
      )
      final_display_rows[[length(final_display_rows) + 1]] <- row_df
    }
  }
}

final_table <- do.call(rbind, final_display_rows)

# ==========================================
# 6. Write to Excel (Formatting)
# ==========================================
# Initialize Workbook
wb <- createWorkbook()
addWorksheet(wb, "Combined Report")

# Write data (NAs will automatically be left blank)
writeData(wb, sheet = "Combined Report", x = final_table, rowNames = FALSE)

# Define Styles
header_style <- createStyle(textDecoration = "bold", border = "Bottom")
assay_header_style <- createStyle(textDecoration = "bold") 
count_style <- createStyle(numFmt = "#,##0") 
wrap_style <- createStyle(wrapText = TRUE)
top_valign_style <- createStyle(valign = "top") # Keeps text top-aligned when rows expand

# Apply Header Style
addStyle(wb, sheet = "Combined Report", style = header_style, rows = 1, cols = 1:ncol(final_table), gridExpand = TRUE)

# Apply bold style to Assay rows
assay_rows <- which(is.na(final_table$`Assay Level`)) + 1 # +1 to account for the header row
for(r in assay_rows) {
  addStyle(wb, sheet = "Combined Report", style = assay_header_style, rows = r, cols = 1)
}

# Apply alignments and wrapping
data_rows <- 2:(nrow(final_table) + 1)
addStyle(wb, sheet = "Combined Report", style = top_valign_style, rows = data_rows, cols = 1:ncol(final_table), gridExpand = TRUE, stack = TRUE)
# Applied wrap_style to columns 1 AND 2 so the Assay Level line breaks work beautifully 
addStyle(wb, sheet = "Combined Report", style = wrap_style, rows = data_rows, cols = 1:2, gridExpand = TRUE, stack = TRUE)

# Apply Count formatting (Cols 3, 4, 5)
addStyle(wb, sheet = "Combined Report", style = count_style, rows = data_rows, cols = 3:5, gridExpand = TRUE, stack = TRUE)

# Adjust Column Widths
setColWidths(wb, sheet = "Combined Report", cols = 1:2, widths = 22) # Ensure width handles the wrapped text for both columns gracefully
setColWidths(wb, sheet = "Combined Report", cols = 3:5, widths = "auto") 

# Save Workbook
saveWorkbook(wb, file.path(out_dir, "agreement_and_confusion_unified.xlsx"), overwrite = TRUE)