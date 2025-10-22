## Given a nano_X.rds, a EM_X.rds, and an output directory,
## return two heatmaps: 1X and 10X.


## Functions ----
sizeMatrix <- function(beta1, beta2, res = 0.01){
  window <- seq(0, 1, res)
  window[1] <- -1
  num_window <- length(window)
  sizeM <- matrix(0, nrow = num_window-1, ncol = num_window-1)
  
  for(i in 1:(num_window-1)){
    print(i)
    flag1 <- (beta1>window[i] & beta1<=window[i+1])
    for(j in 1:(num_window-1)){
      sizeM[i,j] <- sum(flag1 & (beta2>window[j] & beta2<=window[j+1]))
    }
  }
  return(sizeM)
}

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




# Returns a list of two dataframes: nano_sel and em_sel containing only the shared CpG sites
select_data <- function(nano_df, em_df, coverage_cutoff){
  
  nano_df <- nano_df[nano_df$score >= coverage_cutoff, ]
  em_df <- em_df[em_df$score >= coverage_cutoff, ]
  
  em_df$pos <- paste(em_df$chr, em_df$start, sep = "_")
  nano_df$pos <- paste(nano_df$chr, (nano_df$start+1), sep = "_")
  
  share_pos <- em_df$pos[em_df$pos %in% nano_df$pos]
  
  data_em_sel <- em_df[em_df$pos %in% share_pos, ]
  data_nano_sel <- nano_df[nano_df$pos %in% share_pos, ]
  
  if (sum(data_em_sel$pos == data_nano_sel$pos) != nrow(data_em_sel)) {
    stop("sample mismatch")
  }
  
  data_em_sel$beta <- data_em_sel$methy / data_em_sel$score
  data_nano_sel$beta <- (data_nano_sel$NMod + data_nano_sel$NOtherMod) / data_nano_sel$score
  
  return(list(data_nano_sel = data_nano_sel, data_em_sel = data_em_sel))
}





#Retrun a 100 x 3 matrix sizeDF of the number of concordant CpG sites
make_sizeDF <- function(nano_df, em_df, coverage_cutoff) {
  
  data_nano_sel<- nano_df[nano_df$score >= coverage_cutoff, ]
  data_em_sel <- em_df[em_df$score >= coverage_cutoff, ]
  
  
  sizeM <- sizeMatrix(data_em_sel$beta, data_nano_sel$beta, res = 0.01)
  sizeDF <- M2DF(sizeM)
  sizeDF$i <- sizeDF$i / nrow(sizeM)
  sizeDF$j <- sizeDF$j / ncol(sizeM)
  
  return(sizeDF)
}


find_out <- function(x){
  divisors <- 2:x
  smallest <-  divisors[which(x %% divisors == 0)[1]]
  return(smallest + 1)
}



plot_heatmap <- function(sizeDF, nano_sel, em_sel){
  
  gpHeat <- ggplot(sizeDF, aes(x = i, y = j)) +
    geom_raster(aes(fill = log10(size + 1)),  hjust = 0, vjust = 0, interpolate = TRUE) +
    scale_fill_viridis_c(
      option = "D", 
      direction = 1,
      limits = c(0, ceiling(max(log10(sizeDF$size + 1)))), 
      breaks = seq(0, ceiling(max(log10(sizeDF$size + 1))), length.out = find_out(ceiling(max(log10(sizeDF$size + 1))))),
    ) +
    labs(x = "EM-Seq", y = "Nanopore", fill = "log10(Number of Sites)") +
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
  
  
  gpHeat
  
  xDensity <- ggplot(em_sel) +
    geom_density(aes(x = beta*100), fill = "grey", bw = 3) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_pubr() +
    clean_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    rremove("xlab") +
    rremove("ylab")
  
  yDensity <- ggplot(nano_sel) +
    geom_density(aes(x = beta*100), fill = "grey", bw = 3) +
    scale_x_continuous(expand = c(0, 0)) +
    rotate() +
    theme_pubr() +
    clean_theme() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    rremove("xlab") +
    rremove("ylab") +
    scale_y_reverse(expand = c(0, 0))
  
  gp <- ggarrange(NULL, xDensity, yDensity, gpHeat,
                  ncol = 2, nrow = 2, align = "hv",
                  widths = c(1, 6), heights = c(1, 6),
                  common.legend = TRUE, legend = "right"
  )
  
  
  return(gp)
  
}



make_heatmap <- function(nano_df, em_df, coverage_cutoff){
  
  # Select the shared CpG sites between both
  selected <- select_data(nano_df, em_df, coverage_cutoff)
  
  #Make the sizeDF
  size_DF <- make_sizeDF(selected$data_nano_sel, selected$data_em_sel, coverage_cutoff)
  
  #Make the heatmap
  heatmap <- plot_heatmap(size_DF, selected$data_nano_sel, selected$data_em_sel)
  
  return(heatmap)
  
}


heat_filename <- paste0('output/heatmap_s1_10x.png')

png(heat_filename, width = 9, height = 8, units = "in", res = 300)
print(make_heatmap(nano_1, em_1, coverage_cutoff = 10))
dev.off()

heat_filename <- paste0('output/heatmap_s1_1x.png')

png(heat_filename, width = 9, height = 8, units = "in", res = 300)
print(make_heatmap(nano_1, em_1, coverage_cutoff = 1))
dev.off()

### Sample 2

heat_filename <- paste0('output/heatmap_s2_10x.png')

png(heat_filename, width = 9, height = 8, units = "in", res = 300)
print(make_heatmap(nano_2, em_2, coverage_cutoff = 10))
dev.off()

heat_filename <- paste0('output/heatmap_s2_1x.png')

png(heat_filename, width = 9, height = 8, units = "in", res = 300)
print(make_heatmap(nano_2, em_2, coverage_cutoff = 1))
dev.off()

### Sample 3

heat_filename <- paste0('output/heatmap_s3_10x.png')

png(heat_filename, width = 9, height = 8, units = "in", res = 300)
print(make_heatmap(nano_3, em_3, coverage_cutoff = 10))
dev.off()

heat_filename <- paste0('output/heatmap_s3_1x.png')

png(heat_filename, width = 9, height = 8, units = "in", res = 300)
print(make_heatmap(nano_3, em_3, coverage_cutoff = 1))
dev.off()

### Sample 4

heat_filename <- paste0('output/heatmap_s4_10x.png')

png(heat_filename, width = 9, height = 8, units = "in", res = 400)
print(make_heatmap(nano_4, em_4, coverage_cutoff = 10))
dev.off()

heat_filename <- paste0('output/heatmap_s4_1x.png')

png(heat_filename, width = 9, height = 8, units = "in", res = 300)
print(make_heatmap(nano_4, em_4, coverage_cutoff = 1))
dev.off()