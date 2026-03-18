library(tidyverse)
library(GenomicRanges)
library(Biostrings)
library(Rsamtools)
library(patchwork)

# Input:
# {sample} x {threshold} x {region} 
# {s1,s2,s3,s4} x {1X, 10X} x {Genes, transcrtips, ..., whole genome}

# Output:
# Fisher exact test on ONT&EM x capture&missed CpGs for each combo above

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

# Regions
# Load CpG annotation (islands, shores, etc.)
cpg_annotation <- readRDS(file = 'data/reference/cpg_annotation.rds')

# Load gene annotations (hg38 reference)
gene_annotations <- readRDS(file = 'data/reference/gene_annotations.rds')

# Fasta for missing CpGs
fasta_path <- "data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

message("Loading Reference FASTA...")
genome_fasta <- FaFile(fasta_path)

# Function to count unique CpGs in a given GRanges object
# Function to count unique CpGs in a given GRanges object
count_cpgs <- function(gr, fasta) {
  if (is.null(gr)) {
    # Special case: Whole Genome
    whole_genome_seqs <- getSeq(fasta)
    return(sum(vcountPattern("CG", whole_genome_seqs)))
  } else {
    # 1. Merge overlapping regions to prevent double-counting
    gr_reduced <- GenomicRanges::reduce(gr)
    
    # 2. Ensure we only keep chromosomes that exist in the FASTA
    common_chroms <- intersect(seqlevels(gr_reduced), seqlevels(fasta))
    gr_reduced <- keepSeqlevels(gr_reduced, common_chroms, pruning.mode = "coarse")
    
    # 3. Explicitly clamp coordinates instead of using trim()
    # Get valid maximum lengths for the chromosomes currently in gr_reduced
    valid_lengths <- seqlengths(fasta)[as.character(seqnames(gr_reduced))]
    
    # Force start coordinates to be at least 1
    start(gr_reduced) <- pmax(start(gr_reduced), 1L)
    
    # Force end coordinates to not exceed the chromosome's maximum length
    end(gr_reduced) <- pmin(end(gr_reduced), valid_lengths)
    
    # Remove any ranges that became invalid (e.g., width <= 0) after clamping
    gr_reduced <- gr_reduced[width(gr_reduced) > 0]
    
    # 4. Extract sequences securely
    seqs <- getSeq(fasta, gr_reduced)
    
    # 5. Count "CG"
    return(sum(vcountPattern("CG", seqs)))
  }
}


# 1. Define parameters
thresholds <- c(1, 10)
# Create a named list of regions to make looping easy
regions_list <- list(
  "Whole_Genome" = NULL, # Handled as a special case
  "Genes"        = gene_annotations[gene_annotations$type == "gene",],
  "Transcripts"        = gene_annotations[gene_annotations$type == "transcript",],
  "Exons"        = gene_annotations[gene_annotations$type == "exon",],
  "CDS"        = gene_annotations[gene_annotations$type == "CDS",],
  "UTR"        = gene_annotations[gene_annotations$type == "UTR",],
  "Promoter"  = cpg_annotation[cpg_annotation$type == "Promoter", ],
  "CpG Island"  = cpg_annotation[cpg_annotation$type == "CpG_Island", ],
  "Shore"  = cpg_annotation[cpg_annotation$type == "Shore", ],
  "Shelf"  = cpg_annotation[cpg_annotation$type == "Shelf", ],
  "Open Sea"  = cpg_annotation[cpg_annotation$type == "Open_Sea", ]
)

# Pre-calculate and store total CpGs for every region in a list
message("Pre-calculating total CpGs for all regions. This may take a minute...")
total_cpgs_list <- list()
for (region_name in names(regions_list)) {
  total_cpgs_list[[region_name]] <- count_cpgs(regions_list[[region_name]], genome_fasta)
  message(sprintf("  %s: %s total CpGs", region_name, format(total_cpgs_list[[region_name]], big.mark=",")))
}


# Initialize a list to hold the final results
results_list <- list()

# 2. OUTER LOOP: Iterate over samples
for (i in seq_along(file_paths_ONT)) {
  
  # Extract sample ID (e.g., "51", "69")
  sample_id <- str_extract(basename(file_paths_ONT[i]), "\\d+")
  message(sprintf("Loading and processing Sample: %s...", sample_id))
  
  # Load matched data into memory
  ont_data <- readRDS(file_paths_ONT[i])
  em_data  <- readRDS(file_paths_EM[i])
  
  # 3. MIDDLE LOOP: Iterate over thresholds
  for (thresh in thresholds) {
    
    # Filter by coverage 
    ont_filt <- ont_data %>% filter(score >= thresh)
    em_filt  <- em_data %>% filter(score >= thresh) |> mutate(end = start+1)
    
    ont_gr <- makeGRangesFromDataFrame(
      ont_filt, 
      keep.extra.columns = TRUE, 
      seqnames.field = "chr",  
      start.field = "start",   
      end.field = "end"        
    )
    
    em_gr <- makeGRangesFromDataFrame(
      em_filt, 
      keep.extra.columns = TRUE,
      seqnames.field = "chr",
      start.field = "start",
      end.field = "end"
    )
    
    # 4. INNER LOOP: Iterate over regions
    for (region_name in names(regions_list)) {
      
      region_data <- regions_list[[region_name]]
      
      # Step A: Subset your CpG data by the current region
      if (region_name == "Whole_Genome") {
        
        ont_sub_gr <- ont_gr
        em_sub_gr <- em_gr
        
      } else {
        
        # Perform the overlap
        ont_sub_gr <- subsetByOverlaps(ont_gr, region_data)
        em_sub_gr  <- subsetByOverlaps(em_gr, region_data)
        
      }
      
      total_possible_cpgs <- total_cpgs_list[[region_name]]
      
      # Calculate your 2x2 contingency table metrics
      
      ont_unique_gr <- GenomicRanges::reduce(ont_sub_gr)
      em_unique_gr  <- GenomicRanges::reduce(em_sub_gr)
      
      # Calculate counts
      ont_captured <- length(ont_unique_gr)                # ONT captured
      ont_missed <- total_possible_cpgs - ont_captured              # ONT missed
      em_captured <- length(em_unique_gr)              # EM captured
      em_missed <- total_possible_cpgs - em_captured     # EM missed
      
      cont_matrix <- matrix(
        c(em_captured, em_missed,
          ont_captured, ont_missed),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
          "Method" = c("EM_Seq", "Nanopore"),
          "Status" = c("Captured", "Missed")
        )
      )
      
      # Step C: Fisher's Exact Test
      fisher_res <- fisher.test(cont_matrix)
      
      # Step D: Store the results in a tidy format
      res_row <- tibble(
        Sample    = sample_id,
        Threshold = paste0(thresh, "X"),
        Region    = region_name,
        P_Value   = fisher_res$p.value,
        OddsRatio = fisher_res$estimate,
        Total_CpGs = total_possible_cpgs,
        em_captured = em_captured,
        em_missed = em_missed,
        ont_captured = ont_captured,
        ont_missed = ont_missed
      )
      
      results_list[[length(results_list) + 1]] <- res_row
    }
  }
  
  # 5. UNLOAD: Free up memory before loading the next sample
  rm(ont_data, em_data, ont_filt, em_filt)
  gc() # Force garbage collection
}

# 6. Combine all rows into one clean data frame
final_results <- bind_rows(results_list)

# P_value is all zeros...
final_results <- final_results %>%
  mutate(
    # Convert numeric p-values to formatted character strings
    P_Value = case_when(
      P_Value < 0.0001 ~ "< 0.0001",
      TRUE ~ as.character(signif(P_Value, digits = 3)) # Keeps 3 sig figs for larger p-values
    )
  )

write_csv(final_results, file.path(out_dir, "fisher_results_summary.csv"))

final_results <- read.csv("addendum/res/fisher_results_summary.csv")
final_results <- final_results %>%
  mutate(Sample_ID = factor(Sample, 
                            levels = c("51", "69", "168", "266"), 
                            labels = c("A) Sample 1", "B) Sample 2", "C) Sample 3", "D) Sample 4"))) |>
  mutate(Region = factor(Region,
                         levels = c("Genes", "Transcripts", "Exons", "CDS",
                                    "UTR", "CpG Island", "Shore", "Shelf",
                                    "Open Sea", "Promoter", "Whole_Genome")))


# plot

fp <- final_results %>%
  ggplot(aes(x = Region, y = OddsRatio, color = Threshold)) +
  
  # Add a reference line at Odds Ratio = 1 (Equal odds)
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", size = 0.8) +
  
  # Plot points, dodged slightly so 1X and 10X don't overlap
  geom_point(position = position_dodge(width = 0.5), size = 3.5, alpha = 0.9) +
  
  # Facet by Sample
  facet_wrap(~ Sample_ID, ncol = 4) + 
  
  # Use a logarithmic scale for Odds Ratios (best practice)
  scale_y_log10(breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2)) +
  
  # Custom nice colors
  scale_color_manual(values = c("1X" = "#2c7bb6", "10X" = "#d7191c")) +
  
  # LABELS
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Genomic Region",
    # expression() parses the math. ~ adds a space. frac(top, bottom) makes the fraction.
    y = expression("Odds Ratio (EM / ONT) =" ~ frac(EM[det] %*% ONT[miss], EM[miss] %*% ONT[det])), 
    color = "Coverage"
  ) +
  
  # Clean, minimal theme
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    strip.text = element_text(face = "bold", size = 14), # Facet labels
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    legend.position = "bottom"
  )

print(fp)


ggsave(file.path(out_dir, "Suppl_Fig_odds_ratio_EM_Nanopore.png"), plot = fp,
       width = 11, dpi = 300)