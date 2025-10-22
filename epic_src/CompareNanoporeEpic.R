# Compare Nanopore and Epic Array


library(ggplot2)
library(ggpubr)


## Hongyu's email (June 3, 2024)
sample_match <- data.frame(
  nanopore = c("s4059246541", "s4059246577", "s4059249475", "s4059249500"),
  epic = c("RID_4501", "RID_4077", "RID_4466", "RID_4390")
)


## Summary of Epic Array ----
### Load Data
Epic_Manifest <- read.table(
  gzfile("data/EPIC.hg38.manifest.tsv.gz"), 
  header = TRUE, sep = ""
)

Beta_Epic <- read.csv(
  "data/4_individuals_betas 2.csv"
)

cn <- colnames(Beta_Epic)
cn[1] <- "CpG"
colnames(Beta_Epic) <- cn

idx <- stringr::str_detect(Beta_Epic$CpG, pattern = "cg")
Beta_Epic <- Beta_Epic[idx,]
Beta_Epic <- Beta_Epic[Beta_Epic$CpG %in% Epic_Manifest$Probe_ID,]
Beta_Epic <- Beta_Epic[,c("CpG", sample_match$epic)]

idx <- match(Beta_Epic$CpG, Epic_Manifest$Probe_ID)
Epic_Manifest <- Epic_Manifest[idx,]

## 860K CpGs



## Summary of Nanopore Data ----

### Load Data ----
s4059246541 <- readRDS("data/sample_4059246541_analysis/s4059246541_combined.rds")
s4059246577 <- readRDS("data/sample_4059246577_analysis/s4059246577_combined.rds")
s4059249475 <- readRDS("data/sample_4059249475_Analysis/s4059249475_combined.rds")
s4059249500 <- readRDS("data/sample_4059249500_analysis/s4059249500_combined.rds")


mean(s4059246541$score) #22.60951 "A"
mean(s4059249500$score) #28.42207 "B"
mean(s4059246577$score) #29.80659 "C"
mean(s4059249475$score) #25.44275 "D"


mean(
  c(mean(s4059246541$score),
  mean(s4059246577$score),
  mean(s4059249475$score),
  mean(s4059249500$score))
)

#26.57023

write.csv(s4059246541[,c(1,2)], file = "data/CpG_position.csv",
          row.names = FALSE, quote = FALSE)

### Number of CpGs with different coverage  ----
max_coverage <- 30
numCpG_Cov <- matrix(0, nrow = max_coverage, ncol = 4)
for(i in 1:max_coverage){
  numCpG_Cov[i, 1] <- sum(s4059246541$score>=i)
  numCpG_Cov[i, 2] <- sum(s4059246577$score>=i)
  numCpG_Cov[i, 3] <- sum(s4059249475$score>=i)
  numCpG_Cov[i, 4] <- sum(s4059249500$score>=i)
}

numCpG_Cov <- data.frame(
  cbind(1:max_coverage, numCpG_Cov)
)

colnames(numCpG_Cov) <- c("Coverage", sample_match$epic)

write.csv(numCpG_Cov, file = "rlt/tables/numCpG_Cov.csv", row.names = FALSE)

dplot <- reshape2::melt(numCpG_Cov, id.vars = "Coverage", variable.name = "Sample", value.name = "NumCpG")

gp <- ggplot(dplot, aes(x=Coverage, y=NumCpG/10^6, color = Sample)) + 
  geom_point() + 
  #geom_smooth(se = FALSE) + 
  labs(y="Number of CpG (Million)") + 
  scale_color_manual(values = c("#00468BFF", "#42B540FF", "#ED0000FF", "#FDAF91FF")) + 
  scale_y_continuous(limits = c(0, 30)) + 
  theme_pubr() + 
  theme(legend.position = "bottom", legend.title = element_blank())

png("rlt/figures/numCpGCoverage.png", width = 5, height = 4, units = "in", res = 300)
print(gp)
dev.off()


### Beta value distribution ----

#' Figure for Beta value distribution 
#' 
#' @param data_nanopore data frame bedmethyl file
#' @param nSample down sampling CpGs for beta value distribution figure
#' @param file name of file for output
beta_density_all <- function(
    data_nanopore, nSample = NULL, file){
  if(!is.null(nSample)){
    data_nanopore <- data_nanopore[sort(sample(nrow(data_nanopore), nSample)), ]
  }
  
  dplot <- data.frame(
    beta = c(1 - data_nanopore$NCanonical/data_nanopore$score, 
             data_nanopore$NMod/data_nanopore$score,
             data_nanopore$NOtherMod/data_nanopore$score),
    type = factor(c(rep("5mc+5hmc", nrow(data_nanopore)),
                    rep("5mc", nrow(data_nanopore)),
                    rep("5hmc", nrow(data_nanopore))))
  )
  
  gp <- ggplot(dplot) + 
    stat_density(aes(x=beta, color = type), bw = 0.03, geom = "line", position = "identity") + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    labs(x="Beta", y="") + 
    theme_light() +
    theme(legend.position = "bottom", legend.title = element_blank())
  
  png(file, width = 5, height = 4, units = "in", res = 300)
  print(gp)
  dev.off()
}

beta_density_all(s4059246541, nSample = 1000000, file = "rlt/figures/BetaDistAll_RID_4501.png")
beta_density_all(s4059246577, nSample = 1000000, file = "rlt/figures/BetaDistAll_RID_4077.png")
beta_density_all(s4059249475, nSample = 1000000, file = "rlt/figures/BetaDistAll_RID_4466.png")
beta_density_all(s4059249500, nSample = 1000000, file = "rlt/figures/BetaDistAll_RID_4390.png")


# ggplot(data_em_sel) + 
#   geom_density(aes(x=beta), fill = "grey", bw = 3) + 
#   scale_x_continuous(expand = c(0,0)) + 
#   scale_y_continuous(expand = c(0,0)) + 
#   theme_pubr() + 
#   clean_theme() +
#   theme(plot.margin=unit(c(0,0,0,0),"cm")) + 
#   theme(
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank()
#   )  + rremove('xlab') + rremove("ylab")


## Compare Between Nanopore and Epic Array ----

### Original Data ----

Epic_Manifest$chr_pos <- paste(Epic_Manifest$CpG_chrm, Epic_Manifest$CpG_beg, sep = "_")
s4059246541$chr_pos <- paste(s4059246541$chr, s4059246541$start, sep = "_")
s4059246577$chr_pos <- paste(s4059246577$chr, s4059246577$start, sep = "_")
s4059249475$chr_pos <- paste(s4059249475$chr, s4059249475$start, sep = "_")
s4059249500$chr_pos <- paste(s4059249500$chr, s4059249500$start, sep = "_")


#' number of sites with different ranges of beta values. This matrix 
#' will be used for heatmap
#' 
#' @param beta1 beta value 1
#' @param beta2 beta value 2
#' @param res size of the window, there will be 1/res windows 
#' in total
sizeMatrix <- function(beta1, beta2, res = 0.01){
  window <- seq(0, 1, res)
  window[1] <- -1
  num_window <- length(window)
  sizeM <- matrix(0, nrow = num_window-1, ncol = num_window-1)
  
  for(i in 1:(num_window-1)){
    if(i==1){
      flag1 <- (beta1>=0 & beta1<=window[i+1])
    } else {
      flag1 <- (beta1>window[i] & beta1<=window[i+1])
    }
    for(j in 1:(num_window-1)){
      if(j==1){
        flag2 <- beta2 >= 0  & beta2 <= window[j+1]
      } else {
        flag2 <- beta2>window[j] & beta2<=window[j+1]
      }
      sizeM[i,j] <- sum(flag1 & flag2)
    }
  }
  return(sizeM)
}

#' convert matrix to data frame
#' @param M matrix from sizeMatrix
M2DF <- function(M){
  dm <- dim(M)
  size <- dm[1] * dm[2]
  df <- data.frame(
    i = rep(0, size),
    j = rep(0, size),
    size = rep(0, size)
  )
  idx <- 1
  for(i in 1:dm[1]){
    for(j in 1:dm[2]){
      df$i[idx] <- i
      df$j[idx] <- j
      df$size[idx] <- M[i,j]
      idx <- idx + 1
    }
  }
  df
}

min_cov <- 10

#' Correlation between nanopre and epic data
#' 
#' @param data_nanopore nanopore data
#' @param data_epic epic array data
#' @param min_cov minimal coverage for nanopore data
#' @param file_heat name of file to output heatmap
Cor_Nanopore_Epic <- function(
    data_nanopore, data_epic, min_cov = 10,
    file_diff,
    file_heat
){
  data_nanopore <- data_nanopore[data_nanopore$score >= min_cov,]
  num_cpg_nano <- nrow(data_nanopore)
  data_nanopore <- data_nanopore[data_nanopore$chr_pos %in% Epic_Manifest$chr_pos,]
  data_epic <- data_epic[match(data_nanopore$chr_pos, Epic_Manifest$chr_pos),]
  num_overlap <- nrow(data_nanopore)
  colnames(data_epic) <- c("CpG", "beta")
  
  data_nanopore$beta <- 1 - data_nanopore$NCanonical/data_nanopore$score
  
  diff <- data_nanopore$beta - data_epic[,2]
  
  gp <- ggplot(data.frame(diff)) + 
    geom_histogram(
      aes(x=diff), 
      breaks = seq(-1, 1, by = 0.1),
      color="black", fill = "grey") + 
    labs(x="Beta Difference (Nano - Epic)", y = "Count") + 
    theme_pubr()
  
  png(file_diff, width =5, height = 4, units = "in", res = 300)
  print(gp)
  dev.off()
  
  
  
  corr <-  cor(data_nanopore$beta, data_epic[,2], method = "pearson")
  
  sizeM <- sizeMatrix(data_nanopore$beta, data_epic[,2], res = 0.01)
  sizeDF <- M2DF(sizeM)
  sizeDF$i <- sizeDF$i/nrow(sizeM)
  sizeDF$j <- sizeDF$j/ncol(sizeM)
  
  
  # gpHeat <- ggplot(sizeDF, aes(x=i, y=j)) + 
  #   geom_raster(aes(fill = log10(size+1)), hjust = 0,
  #               vjust = 0) + 
  #   #scale_color_gradient(low = "#00468B99", high = "#ED000099") + 
  #   scale_fill_gradient(low = "#00468B", high = "#ED0000", n.breaks = 9)+
  #   labs(x="Nanopore", y="EPIC", fill = "log10(Number of Sites)") + 
  #   #scale_y_continuous(position = "right") + 
  #   scale_x_continuous(limits = c(0,1), expand = c(0,0)) + 
  #   scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  #   #expand_limits(x=c(0, 1), y=c(0, 1)) + 
  #   theme_classic() + 
  #   theme(plot.margin=unit(c(0,0.25,0,0),"cm"), 
  #         legend.title = element_text(angle = -90),
  #         legend.key.height= unit(1, 'cm')) + 
  #   guides(fill = guide_legend(title.position = "right"))
  
  find_out <- function(x){
    divisors <- 2:x
    smallest <-  divisors[which(x %% divisors == 0)[1]]
    return(smallest + 1)
  }
  
  gpHeat <- ggplot(sizeDF, aes(x = i, y = j)) +
    geom_raster(aes(fill = log10(size + 1)),  hjust = 0, vjust = 0, interpolate = TRUE) +
    scale_fill_viridis_c(
      option = "D", 
      direction = 1,
      limits = c(0, ceiling(max(log10(sizeDF$size + 1)))), 
      breaks = seq(0, ceiling(max(log10(sizeDF$size + 1))), length.out = find_out(ceiling(max(log10(sizeDF$size + 1))))),
    ) +
    labs(x = "Nanopore", y = "EPIC", fill = "log10(Number of Sites)") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_fixed() +
    theme_classic() +
    theme(
      plot.margin = unit(c(0, 0.25, 0, 0), "cm"),
      legend.title = element_text(angle = -90, size = 14),
      legend.key.height = unit(1, "cm"),
      axis.title.x = element_text(size = 16),  # Increase X-axis title size
      axis.title.y = element_text(size = 16),  # Increase Y-axis title size
      legend.text = element_text(size = 16)   # Increase legend text size
    ) + 
    guides(fill = guide_colorbar(
      title.position = "right", # Keep the title on the right
      reverse = FALSE, # Flip the legend
      barheight = unit(7, "cm"), # Increase the height of the color bar
      barwidth = unit(0.5, "cm")  # Adjust the width if needed
    ))
  
  
  
  
  
  xDensity <- ggplot(data_nanopore) + 
    geom_density(aes(x=beta), fill = "grey", bw = 0.03) + 
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) + 
    theme_pubr() + 
    clean_theme() +
    theme(plot.margin=unit(c(0,0,0,0),"cm")) + 
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )  + rremove('xlab') + rremove("ylab")
  
  yDensity <- ggplot(data_epic) + 
    geom_density(aes(x=beta), fill = "grey", bw= 0.03)+
    scale_x_continuous(expand = c(0,0)) + 
    rotate() + theme_pubr() + clean_theme()  +
    theme(plot.margin=unit(c(0,0,0,0),"cm")) + 
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) + rremove('xlab') + rremove("ylab") +
    scale_y_reverse(expand = c(0,0))
  
  # xHist <- ggplot(data_nanopore) + 
  #   geom_histogram(aes(x=beta), fill = "grey", breaks = seq(0, 1, by = 0.04)) + 
  #   scale_x_continuous(expand = c(0,0)) + 
  #   scale_y_continuous(expand = c(0,0)) + 
  #   theme_pubr() + 
  #   clean_theme() +
  #   theme(plot.margin=unit(c(0,0,0,0),"cm")) + 
  #   theme(
  #     axis.title.x = element_blank(),
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank()
  #   )  + rremove('xlab') + rremove("ylab")
  # 
  # yHist <- ggplot(data_epic) + 
  #   geom_histogram(aes(x=beta), fill = "grey", breaks = seq(0, 1, by = 0.04)) +
  #   scale_x_continuous(expand = c(0,0)) + 
  #   rotate() + theme_pubr() + clean_theme()  +
  #   theme(plot.margin=unit(c(0,0,0,0),"cm")) + 
  #   theme(
  #     axis.title.y = element_blank(),
  #     axis.text.y = element_blank(),
  #     axis.ticks.y = element_blank()
  #   ) + rremove('xlab') + rremove("ylab") +
  #   scale_y_reverse(expand = c(0,0))
  
  gp <- ggarrange(NULL, xDensity, yDensity, gpHeat, 
                  ncol = 2, nrow = 2,  align = "hv", 
                  widths = c(1, 6), heights = c(1, 6), 
                  common.legend = TRUE, legend = "right")
  
  
  png(file_heat, width = 9, height = 8, units = "in", res = 300)
  print(gp)
  dev.off()
  
  return(c(num_cpg_nano, num_overlap, corr))
}

#### min_cov = 10 ----

cor_beta_RID_4501 <- Cor_Nanopore_Epic(
  s4059246541, Beta_Epic[,c(1,2)], min_cov = 10,
  file_diff = "rlt/figures/Diff_Beta_RID_4501.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4501.png"
)


cor_beta_RID_4077 <- Cor_Nanopore_Epic(
  s4059246577, Beta_Epic[,c(1,3)], min_cov = 10,
  file_diff = "rlt/figures/Diff_Beta_RID_4077.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4077.png"
)


cor_beta_RID_4466 <- Cor_Nanopore_Epic(
  s4059249475, Beta_Epic[,c(1,4)], min_cov = 10,
  file_diff = "rlt/figures/Diff_Beta_RID_4466.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4466.png"
)


cor_beta_RID_4390 <- Cor_Nanopore_Epic(
  s4059249500, Beta_Epic[,c(1,5)], min_cov = 10,
  file_diff = "rlt/figures/Diff_Beta_RID_4390.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4390.png"
)


cor_beta <- data.frame(rbind(
  cor_beta_RID_4501,
  cor_beta_RID_4077,
  cor_beta_RID_4466,
  cor_beta_RID_4390
))

colnames(cor_beta) <- c("num_cpg_nano", "num_overlap", "corr")
cor_beta$pct_cov <- cor_beta$num_overlap/nrow(Epic_Manifest)
write.csv(cor_beta, file = "rlt/tables/cor_beta.csv")


#### min_cov = 1 ----

cor_beta_RID_4501 <- Cor_Nanopore_Epic(
  s4059246541, Beta_Epic[,c(1,2)], min_cov = 1,
  file_diff = "rlt/figures/Diff_Beta_RID_4501_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4501_min_cov1.png"
)


cor_beta_RID_4077 <- Cor_Nanopore_Epic(
  s4059246577, Beta_Epic[,c(1,3)], min_cov = 1,
  file_diff = "rlt/figures/Diff_Beta_RID_4077_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4077_min_cov1.png"
)


cor_beta_RID_4466 <- Cor_Nanopore_Epic(
  s4059249475, Beta_Epic[,c(1,4)], min_cov = 1,
  file_diff = "rlt/figures/Diff_Beta_RID_4466_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4466_min_cov1.png"
)


cor_beta_RID_4390 <- Cor_Nanopore_Epic(
  s4059249500, Beta_Epic[,c(1,5)], min_cov = 1,
  file_diff = "rlt/figures/Diff_Beta_RID_4390_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4390_min_cov1.png"
)


cor_beta <- data.frame(rbind(
  cor_beta_RID_4501,
  cor_beta_RID_4077,
  cor_beta_RID_4466,
  cor_beta_RID_4390
))

colnames(cor_beta) <- c("num_cpg_nano", "num_overlap", "corr")
cor_beta$pct_cov <- cor_beta$num_overlap/nrow(Epic_Manifest)
write.csv(cor_beta, file = "rlt/tables/cor_beta_min_cov1.csv")



### Down Sampling ----

#' Correlation between nanopre and epic data while different level of downsampling nanopore data
#' 
#' @param data_nanopore nanopore data
#' @param data_epic epic array data
#' @param min_cov minimal coverage for nanopore data
#' @param downsampleing_ratio down sampling ratio. default 0.5, down sampling read depth to 50%
#' @parem file_diff name of file to output beta difference histograme
#' @param file_heat name of file to output heatmap
Cor_Nanopore_Epic_downsampling <- function(
    data_nanopore, data_epic, min_cov = 10, downsampleing_ratio = 0.5,
    file_diff, file_heat
){
  data_nanopore$NMod <- rbinom(
    rep(1, nrow(data_nanopore)), 
    data_nanopore$NMod, downsampleing_ratio)
  data_nanopore$NCanonical <- rbinom(
    rep(1, nrow(data_nanopore)), 
    data_nanopore$NCanonical, downsampleing_ratio)
  data_nanopore$NOtherMod <- rbinom(
    rep(1, nrow(data_nanopore)), 
    data_nanopore$NOtherMod, downsampleing_ratio)
  data_nanopore$score <- data_nanopore$NMod + data_nanopore$NCanonical + data_nanopore$NOtherMod
  
  return(Cor_Nanopore_Epic(
    data_nanopore, data_epic, min_cov,file_diff, file_heat
  ))
}

#### min_cov = 10

set.seed(123)

cor_beta_RID_4501 <- Cor_Nanopore_Epic_downsampling(
  s4059246541, Beta_Epic[,c(1,2)], min_cov = 10, downsampleing_ratio = 15/23.25,
  file_diff = "rlt/figures/Diff_Beta_RID_4501_down15.png", 
  file_heat = "rlt/figures/Cor_Beta_RID_4501_down15.png"
)


cor_beta_RID_4077 <- Cor_Nanopore_Epic_downsampling(
  s4059246577, Beta_Epic[,c(1,3)], min_cov = 10, downsampleing_ratio = 15/29.32,
  file_diff = "rlt/figures/Diff_Beta_RID_4077_down15.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4077_down15.png"
)


cor_beta_RID_4466 <- Cor_Nanopore_Epic_downsampling(
  s4059249475, Beta_Epic[,c(1,4)], min_cov = 10, downsampleing_ratio = 15/25.43,
  file_diff = "rlt/figures/Diff_Beta_RID_4466_down15.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4466_down15.png"
)


cor_beta_RID_4390 <- Cor_Nanopore_Epic_downsampling(
  s4059249500, Beta_Epic[,c(1,5)], min_cov = 10, downsampleing_ratio = 15/27.35,
  file_diff = "rlt/figures/Diff_Beta_RID_4390_down15.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4390_down15.png"
)


cor_beta <- data.frame(rbind(
  cor_beta_RID_4501,
  cor_beta_RID_4077,
  cor_beta_RID_4466,
  cor_beta_RID_4390
))

colnames(cor_beta) <- c("num_cpg_nano", "num_overlap", "corr")
cor_beta$pct_cov <- cor_beta$num_overlap/nrow(Epic_Manifest)
write.csv(cor_beta, file = "rlt/tables/cor_beta_down15.csv")




cor_beta_RID_4501 <- Cor_Nanopore_Epic_downsampling(
  s4059246541, Beta_Epic[,c(1,2)], min_cov = 10, downsampleing_ratio = 12.5/23.25,
  file_diff = "rlt/figures/Diff_Beta_RID_4501_down12.5.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4501_down12.5.png"
)


cor_beta_RID_4077 <- Cor_Nanopore_Epic_downsampling(
  s4059246577, Beta_Epic[,c(1,3)], min_cov = 10, downsampleing_ratio = 12.5/29.32,
  file_diff = "rlt/figures/Diff_Beta_RID_4077_down12.5.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4077_down12.5.png"
)


cor_beta_RID_4466 <- Cor_Nanopore_Epic_downsampling(
  s4059249475, Beta_Epic[,c(1,4)], min_cov = 10, downsampleing_ratio = 12.5/25.43,
  file_diff = "rlt/figures/Diff_Beta_RID_4466_down12.5.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4466_down12.5.png"
)


cor_beta_RID_4390 <- Cor_Nanopore_Epic_downsampling(
  s4059249500, Beta_Epic[,c(1,5)], min_cov = 10, downsampleing_ratio = 12.5/27.35,
  file_diff = "rlt/figures/Diff_Beta_RID_4390_down12.5.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4390_down12.5.png"
)


cor_beta <- data.frame(rbind(
  cor_beta_RID_4501,
  cor_beta_RID_4077,
  cor_beta_RID_4466,
  cor_beta_RID_4390
))

colnames(cor_beta) <- c("num_cpg_nano", "num_overlap", "corr")
cor_beta$pct_cov <- cor_beta$num_overlap/nrow(Epic_Manifest)
write.csv(cor_beta, file = "rlt/tables/cor_beta_down12.5.csv")



cor_beta_RID_4501 <- Cor_Nanopore_Epic_downsampling(
  s4059246541, Beta_Epic[,c(1,2)], min_cov = 10, downsampleing_ratio = 10/23.25,
  file_diff = "rlt/figures/Diff_Beta_RID_4501_down10.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4501_down10.png"
)


cor_beta_RID_4077 <- Cor_Nanopore_Epic_downsampling(
  s4059246577, Beta_Epic[,c(1,3)], min_cov = 10, downsampleing_ratio = 10/29.32,
  file_diff = "rlt/figures/Diff_Beta_RID_4077_down10.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4077_down10.png"
)


cor_beta_RID_4466 <- Cor_Nanopore_Epic_downsampling(
  s4059249475, Beta_Epic[,c(1,4)], min_cov = 10, downsampleing_ratio = 10/25.43,
  file_diff = "rlt/figures/Diff_Beta_RID_4466_down10.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4466_down10.png"
)


cor_beta_RID_4390 <- Cor_Nanopore_Epic_downsampling(
  s4059249500, Beta_Epic[,c(1,5)], min_cov = 10, downsampleing_ratio = 10/27.35,
  file_diff = "rlt/figures/Diff_Beta_RID_4390_down10.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4390_down10.png"
)


cor_beta <- data.frame(rbind(
  cor_beta_RID_4501,
  cor_beta_RID_4077,
  cor_beta_RID_4466,
  cor_beta_RID_4390
))

colnames(cor_beta) <- c("num_cpg_nano", "num_overlap", "corr")
cor_beta$pct_cov <- cor_beta$num_overlap/nrow(Epic_Manifest)
write.csv(cor_beta, file = "rlt/tables/cor_beta_down10.csv")




#### min_cov = 1 ----

set.seed(123)

cor_beta_RID_4501 <- Cor_Nanopore_Epic_downsampling(
  s4059246541, Beta_Epic[,c(1,2)], min_cov = 1, downsampleing_ratio = 15/23.25,
  file_diff = "rlt/figures/Diff_Beta_RID_4501_down15_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4501_down15_min_cov1.png"
)


cor_beta_RID_4077 <- Cor_Nanopore_Epic_downsampling(
  s4059246577, Beta_Epic[,c(1,3)], min_cov = 1, downsampleing_ratio = 15/29.32,
  file_diff = "rlt/figures/Diff_Beta_RID_4077_down15_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4077_down15_min_cov1.png"
)


cor_beta_RID_4466 <- Cor_Nanopore_Epic_downsampling(
  s4059249475, Beta_Epic[,c(1,4)], min_cov = 1, downsampleing_ratio = 15/25.43,
  file_diff = "rlt/figures/Diff_Beta_RID_4466_down15_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4466_down15_min_cov1.png"
)


cor_beta_RID_4390 <- Cor_Nanopore_Epic_downsampling(
  s4059249500, Beta_Epic[,c(1,5)], min_cov = 1, downsampleing_ratio = 15/27.35,
  file_diff = "rlt/figures/Diff_Beta_RID_4390_down15_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4390_down15_min_cov1.png"
)


cor_beta <- data.frame(rbind(
  cor_beta_RID_4501,
  cor_beta_RID_4077,
  cor_beta_RID_4466,
  cor_beta_RID_4390
))

colnames(cor_beta) <- c("num_cpg_nano", "num_overlap", "corr")
cor_beta$pct_cov <- cor_beta$num_overlap/nrow(Epic_Manifest)
write.csv(cor_beta, file = "rlt/tables/cor_beta_down15_min_cov1.csv")




cor_beta_RID_4501 <- Cor_Nanopore_Epic_downsampling(
  s4059246541, Beta_Epic[,c(1,2)], min_cov = 1, downsampleing_ratio = 12.5/23.25,
  file_diff = "rlt/figures/Diff_Beta_RID_4501_down12.5_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4501_down12.5_min_cov1.png"
)


cor_beta_RID_4077 <- Cor_Nanopore_Epic_downsampling(
  s4059246577, Beta_Epic[,c(1,3)], min_cov = 1, downsampleing_ratio = 12.5/29.32,
  file_diff = "rlt/figures/Diff_Beta_RID_4077_down12.5_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4077_down12.5_min_cov1.png"
)


cor_beta_RID_4466 <- Cor_Nanopore_Epic_downsampling(
  s4059249475, Beta_Epic[,c(1,4)], min_cov = 1, downsampleing_ratio = 12.5/25.43,
  file_diff = "rlt/figures/Diff_Beta_RID_4466_down12.5_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4466_down12.5_min_cov1.png"
)


cor_beta_RID_4390 <- Cor_Nanopore_Epic_downsampling(
  s4059249500, Beta_Epic[,c(1,5)], min_cov = 1, downsampleing_ratio = 12.5/27.35,
  file_diff = "rlt/figures/Diff_Beta_RID_4390_down12.5_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4390_down12.5_min_cov1.png"
)


cor_beta <- data.frame(rbind(
  cor_beta_RID_4501,
  cor_beta_RID_4077,
  cor_beta_RID_4466,
  cor_beta_RID_4390
))

colnames(cor_beta) <- c("num_cpg_nano", "num_overlap", "corr")
cor_beta$pct_cov <- cor_beta$num_overlap/nrow(Epic_Manifest)
write.csv(cor_beta, file = "rlt/tables/cor_beta_down12.5_min_cov1.csv")



cor_beta_RID_4501 <- Cor_Nanopore_Epic_downsampling(
  s4059246541, Beta_Epic[,c(1,2)], min_cov = 1, downsampleing_ratio = 10/23.25,
  file_diff = "rlt/figures/Diff_Beta_RID_4501_down10_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4501_down10_min_cov1.png"
)


cor_beta_RID_4077 <- Cor_Nanopore_Epic_downsampling(
  s4059246577, Beta_Epic[,c(1,3)], min_cov = 1, downsampleing_ratio = 10/29.32,
  file_diff = "rlt/figures/Diff_Beta_RID_4077_down10_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4077_down10_min_cov1.png"
)


cor_beta_RID_4466 <- Cor_Nanopore_Epic_downsampling(
  s4059249475, Beta_Epic[,c(1,4)], min_cov = 1, downsampleing_ratio = 10/25.43,
  file_diff = "rlt/figures/Diff_Beta_RID_4466_down10_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4466_down10_min_cov1.png"
)


cor_beta_RID_4390 <- Cor_Nanopore_Epic_downsampling(
  s4059249500, Beta_Epic[,c(1,5)], min_cov = 1, downsampleing_ratio = 10/27.35,
  file_diff = "rlt/figures/Diff_Beta_RID_4390_down10_min_cov1.png",
  file_heat = "rlt/figures/Cor_Beta_RID_4390_down10_min_cov1.png"
)


cor_beta <- data.frame(rbind(
  cor_beta_RID_4501,
  cor_beta_RID_4077,
  cor_beta_RID_4466,
  cor_beta_RID_4390
))

colnames(cor_beta) <- c("num_cpg_nano", "num_overlap", "corr")
cor_beta$pct_cov <- cor_beta$num_overlap/nrow(Epic_Manifest)
write.csv(cor_beta, file = "rlt/tables/cor_beta_down10_min_cov1.csv")

