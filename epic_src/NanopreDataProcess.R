
# functions ------


#' combine hap1, hap2, and ungrouped data into single bedmethyl file with 5mc and 5hmc combined
#'
#' @param bedmethyl_file_hap1 bedmethyl file from haplotype 1
#' @param bedmethyl_file_hap2 bedmethyl file from haplotype 2
#' @param bedmethyl_file_ungrouped ungrouped bedmethyl file
#' @param filename_out name of the file to output results. If it is NULL, return the resutls. 
mergeHapsAllMethylation <- function(
    bedmethyl_file_hap1, bedmethyl_file_hap2, 
    bedmethyl_file_ungrouped, filename_out = NULL
){
  cn <- c(
    "chr", "start", "end", "modiCode", "score", "strand", "start2", "end2",
    "color", "NVCov", "percentModi", "NMod", "NCanonical", "NOtherMod", 
    "NDelete", "NFail", "NDiff", "NNoCall"
  )
  
  cn_sel <- c(
    "chr", "start", "end", "modiCode", "score", "strand", "percentModi", 
    "NMod", "NCanonical", "NOtherMod"
  )
  
  bedmethy_hap1 <- read.table(
    gzfile(bedmethyl_file_hap1), 
    header = FALSE, sep = "", col.names = cn
  )   
  
  bedmethy_hap2 <- read.table(
    gzfile(bedmethyl_file_hap2),
    header = FALSE, sep = "", col.names = cn
  )
  
  bedmethy_ungrouped <- read.table(
    gzfile(bedmethyl_file_ungrouped),
    header = FALSE, sep = "", col.names = cn
  )
  
  bedmethy_hap1 <- bedmethy_hap1[,cn_sel]
  bedmethy_hap2 <- bedmethy_hap2[,cn_sel]
  bedmethy_ungrouped <- bedmethy_ungrouped[,cn_sel]
  
  # print(table(bedmethy_hap1$modiCode))
  # print(table(bedmethy_hap2$modiCode))
  # print(table(bedmethy_ungrouped$modiCode))
  
  # not including "MT" here, CpGs on MT are all in ungrouped file
  chr_all <- paste0("chr", c(1:22, "X", "Y"))
  rlt_final <- NULL
  for(chr in chr_all){
    print(chr)
    h1 <- bedmethy_hap1[bedmethy_hap1$chr == chr & bedmethy_hap1$modiCode == "m", ]
    h2 <- bedmethy_hap2[bedmethy_hap2$chr == chr & bedmethy_hap2$modiCode == "m", ]
    ug <- bedmethy_ungrouped[bedmethy_ungrouped$chr == chr & bedmethy_ungrouped$modiCode == "m", ]
    
    start_all <- sort(unique(c(h1$start, h2$start, ug$start)))
    
    rlt <- data.frame(
      chr = chr, 
      start = start_all,
      end = 0,
      modeCode = "m",
      score = 0,
      NMod = 0,
      NCanonical = 0,
      NOtherMod = 0,
      stringsAsFactors = FALSE
    )
    
    # add h1
    idx <- match(h1$start, rlt$start)
    rlt$end[idx] <- h1$end
    rlt$score[idx] <- rlt$score[idx] + h1$score
    rlt$NMod[idx] <- rlt$NMod[idx] + h1$NMod
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + h1$NCanonical
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + h1$NOtherMod
    
    # add h2
    idx <- match(h2$start, rlt$start)
    rlt$end[idx] <- h2$end
    rlt$score[idx] <- rlt$score[idx] + h2$score
    rlt$NMod[idx] <- rlt$NMod[idx] + h2$NMod
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + h2$NCanonical
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + h2$NOtherMod
    
    # add ug
    idx <- match(ug$start, rlt$start)
    rlt$end[idx] <- ug$end
    rlt$score[idx] <- rlt$score[idx] + ug$score
    rlt$NMod[idx] <- rlt$NMod[idx] + ug$NMod
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + ug$NCanonical
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + ug$NOtherMod
    
    if(chr == "chr1"){
      rlt_final <- rlt
    } else {
      rlt_final <- rbind(rlt_final, rlt)
    }
  }
  
  if(is.null(filename_out)){
    return(rlt_final)
  }
  
  saveRDS(rlt_final, file = filename_out)
}



mergeHapsAllMethylation(
  "data/sample_4059246541_analysis/s4059246541_1.wf_mods.bedmethyl.gz",
  "data/sample_4059246541_analysis/s4059246541_2.wf_mods.bedmethyl.gz",
  "data/sample_4059246541_analysis/s4059246541_ungrouped.wf_mods.bedmethyl.gz",
  "data/sample_4059246541_analysis/s4059246541_combined.rds"
)

mergeHapsAllMethylation(
  "data/sample_4059246577_analysis/s4059246577_1.wf_mods.bedmethyl.gz",
  "data/sample_4059246577_analysis/s4059246577_2.wf_mods.bedmethyl.gz",
  "data/sample_4059246577_analysis/s4059246577_ungrouped.wf_mods.bedmethyl.gz",
  "data/sample_4059246577_analysis/s4059246577_combined.rds"
)

mergeHapsAllMethylation(
  "data/sample_4059249475_Analysis/s4059249475_1.wf_mods.bedmethyl.gz",
  "data/sample_4059249475_Analysis/s4059249475_2.wf_mods.bedmethyl.gz",
  "data/sample_4059249475_Analysis/s4059249475_ungrouped.wf_mods.bedmethyl.gz",
  "data/sample_4059249475_Analysis/s4059249475_combined.rds"
)

mergeHapsAllMethylation(
  "data/sample_4059249500_analysis/s4059249500_1.wf_mods.bedmethyl.gz",
  "data/sample_4059249500_analysis/s4059249500_2.wf_mods.bedmethyl.gz",
  "data/sample_4059249500_analysis/s4059249500_ungrouped.wf_mods.bedmethyl.gz",
  "data/sample_4059249500_analysis/s4059249500_combined.rds"
)




# 6mA -----

bedmethyl_file_hap1 <- "data/sample_4059246541_analysis/s4059246541_1.wf_mods.bedmethyl.gz"
bedmethyl_file_hap2 <- "data/sample_4059246541_analysis/s4059246541_2.wf_mods.bedmethyl.gz"
bedmethyl_file_ungrouped <- "data/sample_4059246541_analysis/s4059246541_ungrouped.wf_mods.bedmethyl.gz"
filename_out <- "data/sample_4059246541_analysis/s4059246541_combined_6mA.rds"


mergeHapsAllMethylation <- function(
    bedmethyl_file_hap1, bedmethyl_file_hap2, 
    bedmethyl_file_ungrouped, filename_out = NULL
){
  cn <- c(
    "chr", "start", "end", "modiCode", "score", "strand", "start2", "end2",
    "color", "NVCov", "percentModi", "NMod", "NCanonical", "NOtherMod", 
    "NDelete", "NFail", "NDiff", "NNoCall"
  )
  
  cn_sel <- c(
    "chr", "start", "end", "modiCode", "score", "strand", "percentModi", 
    "NMod", "NCanonical", "NOtherMod"
  )
  
  bedmethy_hap1 <- read.table(
    gzfile(bedmethyl_file_hap1), 
    header = FALSE, sep = "", col.names = cn
  )   
  
  bedmethy_hap2 <- read.table(
    gzfile(bedmethyl_file_hap2),
    header = FALSE, sep = "", col.names = cn
  )
  
  bedmethy_ungrouped <- read.table(
    gzfile(bedmethyl_file_ungrouped),
    header = FALSE, sep = "", col.names = cn
  )
  
  bedmethy_hap1 <- bedmethy_hap1[,cn_sel]
  bedmethy_hap2 <- bedmethy_hap2[,cn_sel]
  bedmethy_ungrouped <- bedmethy_ungrouped[,cn_sel]
  
  # print(table(bedmethy_hap1$modiCode))
  # print(table(bedmethy_hap2$modiCode))
  # print(table(bedmethy_ungrouped$modiCode))
  
  # not including "MT" here, CpGs on MT are all in ungrouped file
  chr_all <- paste0("chr", c(1:22, "X", "Y"))
  rlt_final <- NULL
  for(chr in chr_all){
    print(chr)
    h1 <- bedmethy_hap1[bedmethy_hap1$chr == chr & bedmethy_hap1$modiCode == "a", ]
    h2 <- bedmethy_hap2[bedmethy_hap2$chr == chr & bedmethy_hap2$modiCode == "a", ]
    ug <- bedmethy_ungrouped[bedmethy_ungrouped$chr == chr & bedmethy_ungrouped$modiCode == "a", ]
    
    start_all <- sort(unique(c(h1$start, h2$start, ug$start)))
    
    rlt <- data.frame(
      chr = chr, 
      start = start_all,
      end = 0,
      modeCode = "a",
      score = 0,
      NMod = 0,
      NCanonical = 0,
      NOtherMod = 0,
      stringsAsFactors = FALSE
    )
    
    # add h1
    idx <- match(h1$start, rlt$start)
    rlt$end[idx] <- h1$end
    rlt$score[idx] <- rlt$score[idx] + h1$score
    rlt$NMod[idx] <- rlt$NMod[idx] + h1$NMod
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + h1$NCanonical
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + h1$NOtherMod
    
    # add h2
    idx <- match(h2$start, rlt$start)
    rlt$end[idx] <- h2$end
    rlt$score[idx] <- rlt$score[idx] + h2$score
    rlt$NMod[idx] <- rlt$NMod[idx] + h2$NMod
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + h2$NCanonical
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + h2$NOtherMod
    
    # add ug
    idx <- match(ug$start, rlt$start)
    rlt$end[idx] <- ug$end
    rlt$score[idx] <- rlt$score[idx] + ug$score
    rlt$NMod[idx] <- rlt$NMod[idx] + ug$NMod
    rlt$NCanonical[idx] <- rlt$NCanonical[idx] + ug$NCanonical
    rlt$NOtherMod[idx] <- rlt$NOtherMod[idx] + ug$NOtherMod
    
    if(chr == "chr1"){
      rlt_final <- rlt
    } else {
      rlt_final <- rbind(rlt_final, rlt)
    }
  }
  
  if(is.null(filename_out)){
    return(rlt_final)
  }
  
  saveRDS(rlt_final, file = filename_out)
}


