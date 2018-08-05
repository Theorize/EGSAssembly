library(ggplot2)
library(scales)
library(reshape2)
library(ggthemes)
library(plyr)
library(wesanderson)
setwd("/Users/isla/Documents/University/Masters/Thesis/Code/truncated-suffix-array/DataResults/")



plot_data_dist <- function(data_file, plot_file, do_legend=FALSE) {
  base_dist <- read.table(data_file, header = TRUE, sep = ",")
  base_dist$Outcome <- ifelse(base_dist$base == ".", "Ambiguous",
                              ifelse(base_dist$base == "n", "Too little data", "Accepted"))
  sorted = apply(base_dist[, 1:4], 1, sort)
  base_dist$max <- sorted[4,]
  base_dist$Major <- sorted[3,]
  base_dist$Minor <- sorted[2,]
  
  outcomes_prop <- table(base_dist$Outcome)
  outcomes_prop <- data.frame(round(prop.table(outcomes_prop), digits=7)*100)
  print(paste("Major is 0:", (nrow(subset(base_dist, Major == 0))/nrow(base_dist))*100,"%"))
  print(paste("Accepted (", outcomes_prop[2][[1]][1], "%)", sep=""))
  print(paste("Ambiguous (", outcomes_prop[2][[1]][2], "%)", sep=""))
  print(paste("Not Enough Data (", outcomes_prop[2][[1]][3], "%)", sep=""))
  
  
  accepted_label <- "Accepted"
  ambiguous_label <- "Ambiguous"
  small_label <- "Suffix Interval Under-Populated"
  
  if(folder == "human30") {
    m_limit = 0.25
  } else if (folder == "human60") {
    m_limit = 0.25
  } else {
    m_limit = 0.1
  }
  
  abc <- ggplot(base_dist, aes(x=Major, y=Minor, color=Outcome)) +
         geom_point(alpha=1/10) +  xlim(0, 0.5) + ylim(0, 0.4) + coord_equal(ratio=1) +
         geom_vline(aes(xintercept = m_limit),
                     color="#0B3C49", linetype="dotted", size=0.5) +
         geom_hline(aes(yintercept = m_limit),
               color="#0B3C49", linetype="dotted", size=0.5) +
         scale_color_manual(labels = c("Accepted" = accepted_label, "Ambiguous" = ambiguous_label,
                                       "Too little data" = small_label),
                            values=c("Accepted" = "#0A369D", "Ambiguous" = "#731963",
                                     "Too little data" = "#BBBE64")) +
        theme_tufte()
  
  if (do_legend) {
    abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"),
                       panel.grid.major=element_line(size = 0.1, colour = "gray"),
                       legend.position="bottom", legend.title=element_blank()) +
      guides(colour = guide_legend(override.aes = list(alpha = 1)))
    ggsave(filename = plot_file, abc, width = 11.5, height = 12.5, units = "cm")
  } else {
    abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"),
                       panel.grid.major=element_line(size = 0.1, colour = "gray"),
                       legend.position="none")
    ggsave(filename = plot_file, abc, width = 11.5, height = 11.5, units = "cm")
  }
      
}








plot_data_dist_with_overlay <- function(data_file, plot_file, do_legend=FALSE) {
  base_dist <- read.table(data_file, header = TRUE, sep = ",")
  base_dist$Outcome <- ifelse(base_dist$base == ".", "Ambiguous",
                              ifelse(base_dist$base == "n", "Too little data", "Accepted"))
  sorted = apply(base_dist[, 1:4], 1, sort)
  base_dist$max <- sorted[4,]
  base_dist$Major <- sorted[3,]
  base_dist$Minor <- sorted[2,]
  
  outcomes_prop <- table(base_dist$Outcome)
  outcomes_prop <- data.frame(round(prop.table(outcomes_prop), digits=7)*100)
  print(paste("Major is 0:", (nrow(subset(base_dist, Major == 0))/nrow(base_dist))*100,"%"))
  print(paste("Accepted (", outcomes_prop[2][[1]][1], "%)", sep=""))
  print(paste("Ambiguous (", outcomes_prop[2][[1]][2], "%)", sep=""))
  print(paste("Not Enough Data (", outcomes_prop[2][[1]][3], "%)", sep=""))
  
  
  accepted_label <- "Accepted"
  ambiguous_label <- "Ambiguous"
  small_label <- "Suffix Interval Under-Populated"
  
  if(folder == "human30") {
    m_limit = 0.25
  } else if (folder == "human60") {
    m_limit = 0.25
  } else {
    m_limit = 0.1
  }
  
  abc <- ggplot(base_dist, aes(x=Major, y=Minor, color=Outcome)) +
    geom_point(alpha=1/10) + coord_equal(ratio=1) +
    scale_color_manual(labels = c("Accepted" = accepted_label, "Ambiguous" = ambiguous_label,
                                  "Too little data" = small_label),
                       values=c("Accepted" = "#0A369D", "Ambiguous" = "#731963",
                                "Too little data" = "#BBBE64")) +
    scale_x_continuous(breaks=c(0, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 2/5, 2/7, 3/7, 3/8),
                       labels = c("0", expression(frac(1, 8)), expression(frac(1, 7)), expression(frac(1, 6)), expression(frac(1, 5)), expression(frac(1, 4)),
                                  expression(frac(2, 7)), expression(frac(1, 3)), expression(frac(3, 8)), expression(frac(2, 5)), expression(frac(1, 2))),
                      limits = c(0, 0.5), minor_breaks = NULL) +
     scale_y_continuous(breaks=c(0, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 2/7, 3/7),
                        labels = NULL,
                        limits = c(0, 0), minor_breaks = NULL) + 
    theme_tufte() + theme(axis.line = element_line(size = 0.25, colour = "black"),
                          axis.text.x = element_text(size=8),
                       panel.grid.major=element_line(size = 0.125, colour = "black"),
                       legend.position="none")
    ggsave(filename = plot_file, abc, width = 11.5, height = 12.5, units = "cm")
    
    abc
}

folder = "human60"
kmersize = 15
data_file <- paste("./", folder , "/results_for_", kmersize, "/base_dist.txt", sep="")
plot_file <- file.path(paste("../../../WriteUp/plots/", folder, "/base_dist_with_grid_", kmersize, ".png", sep=""))
plot_data_dist_with_overlay(data_file, plot_file)



  
  
  
  
  
create_contig_len_hist <- function(binsize, folder, kmersize) {
  # Collect data.
  data_file <- paste("./", folder , "/results_for_", kmersize, "/contig_length.txt", sep="")
  contig_len <- read.table(data_file, header = FALSE, sep = ",", col.names=c("Length"))
  contig_len <- subset(contig_len, Length > kmersize + 1)
  
  # Create the plot
  hist_plot <- ggplot(contig_len, aes(x = Length), mai=c(0.25, 0.25, 0.25, 0.25)) +
               geom_histogram(binwidth=binsize, fill = "#0B3C49") +
    xlab(paste("Contig Length (Binsize: ", binsize, ")", sep="")) + xlim(0, NA) +
    ylab("Count (Log1p Scale)") +
    theme_tufte() +
    theme(axis.line = element_line(size = 0.25, colour = "black"),
          panel.grid.major=element_line(size = 0.1, colour = "gray"),
          panel.grid.minor=element_line(size = 0.05, colour = "gray")) +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000),
                       minor_breaks = c(2,3,4,6,7,8,9,20,30,40,60,70,80,90,200,300,400,600,700,800,
                                        900,2000,3000,4000,6000,7000,8000,9000),
                       trans="log1p")
  
  # Save the plot
  plot_file <- file.path(paste("../../../WriteUp/plots/", folder, "/contig_length_hist", kmersize, ".png", sep=""))
  ggsave(filename = plot_file, hist_plot, width = 11.5, height = 11.5, units = "cm")
}







create_contig_len_hist_human <- function(binsize, folder, kmersize) {
  # Collect data.
  data_file <- paste("./", folder , "/results_for_", kmersize, "/contig_length2.txt", sep="")
  contig_len <- read.table(data_file, header = FALSE, sep = ",", col.names=c("Length", "ER"))
  data_file <- paste("./", folder , "/results_for_", kmersize, "/matches.txt", sep="")
  matches <- read.table(data_file, header = FALSE, sep = ";", col.names=c("Repeat"))
  contig_len$Repeat = matches$Repeat
  
  
  repeat_no_check <- subset(contig_len, Length > kmersize*3)
  repeat_no_check$Check = ifelse(repeat_no_check$ER == repeat_no_check$Repeat, 0,
                                 ifelse(repeat_no_check$Repeat == -1, 0, 1))
  print(paste("Consistency check: ", sum(repeat_no_check$Check)/nrow(repeat_no_check)))
  
  contig_len$Prevalence <- ifelse(contig_len$Repeat == 1, "Non-Repetitive", 
                                  ifelse(contig_len$Repeat <= 4, "Repetitive (<5)",
                                         paste("Repetitive (<", max(contig_len$Repeat)+1, ")", sep="")))
  contig_len <- subset(contig_len, Length > kmersize*3)
  
  print(paste(folder, kmersize, "min:", min(contig_len$Repeat), "max:", max(contig_len$Repeat)))
  
  # Create the plot
  hist_plot <- ggplot(contig_len, aes(x = Length)) + geom_histogram(binwidth=binsize, fill = "#0B3C49") +
    xlab(paste("Contig Length (Binsize: ", binsize, ")", sep="")) + xlim(0, NA) +
    ylab("Count (Log1p Scale)") +
    theme_tufte() +
    theme(axis.line = element_line(size = 0.25, colour = "black"),
          panel.grid.major=element_line(size = 0.1, colour = "gray"),
          panel.grid.minor=element_line(size = 0.05, colour = "gray")) +
    scale_y_continuous(breaks=c(0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000),
                       minor_breaks = c(2,3,4,6,7,8,9,20,30,40,60,70,80,90,200,300,400,600,700,800,
                                        900,2000,3000,4000,6000,7000,8000,9000),
                       trans="log1p")
  # Save the plot
  plot_file <- file.path(paste("../../../WriteUp/plots/", folder, "/contig_length_hist", kmersize, ".png", sep=""))
  ggsave(filename = plot_file, hist_plot, width = 11.5, height = 11.5, units = "cm")
  
  # For each prevalence, create plots.
  prev_non <- subset(contig_len, Repeat == 1)
  if (nrow(prev_non)>0) {
    hist_plot <- ggplot(prev_non, aes(x = Length)) + geom_histogram(binwidth=binsize, fill = "#0B3C49") +
      xlab(paste("Contig Length (Binsize: ", binsize, ")", sep="")) + xlim(0, NA) +
      ylab("Count (Log1p Scale)") +
      theme_tufte() +
      theme(axis.line = element_line(size = 0.25, colour = "black"),
            panel.grid.major=element_line(size = 0.1, colour = "gray"),
            panel.grid.minor=element_line(size = 0.05, colour = "gray")) +
      scale_y_continuous(breaks=c(0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000),
                         minor_breaks = c(2,3,4,6,7,8,9,20,30,40,60,70,80,90,200,300,400,600,700,800,
                                          900,2000,3000,4000,6000,7000,8000,9000),
                         trans="log1p")
    plot_file <- file.path(paste("../../../WriteUp/plots/", folder,
                                 "/non_repetitive_contig_length_hist", kmersize, ".png", sep=""))
    ggsave(filename = plot_file, hist_plot, width = 11.5, height = 11.5, units = "cm")
  }
  
  prev_rep <- subset(contig_len, Repeat > 1)
  if (nrow(prev_rep)>0){
    if (kmersize < 36) { reduced_binsize = binsize/5 } else { reduced_binsize = binsize/10 }
    hist_plot <- ggplot(prev_rep, aes(x = Length)) + geom_histogram(binwidth=reduced_binsize, fill = "#0B3C49") +
      xlab(paste("Contig Length (Binsize: ", reduced_binsize, ")", sep="")) + xlim(0, NA) +
      ylab("Count (Log1p Scale)") +
      theme_tufte() +
      theme(axis.line = element_line(size = 0.25, colour = "black"),
            panel.grid.major=element_line(size = 0.1, colour = "gray"),
            panel.grid.minor=element_line(size = 0.05, colour = "gray")) +
      scale_y_continuous(breaks=c(0, 1, 5, 10, 50, 100, 500, 1000, 5000, 10000),
                         minor_breaks = c(2,3,4,6,7,8,9,20,30,40,60,70,80,90,200,300,400,600,700,800,
                                          900,2000,3000,4000,6000,7000,8000,9000),
                         trans="log1p")
    plot_file <- file.path(paste("../../../WriteUp/plots/", folder,
                                 "/repetitive_contig_length_hist", kmersize, ".png", sep=""))
    ggsave(filename = plot_file, hist_plot, width = 11.5, height = 11.5, units = "cm")
  }
}








###### Contig End Reasons

plot_contig_end_reasons <- function(folder, reverse) {
  kmersize = 15
  if (reverse) {
    data_file <- paste("./", folder , "/results_for_15/reverse_contig_end_reasons.txt", sep="")
  } else {
    data_file <- paste("./", folder , "/results_for_15/contig_end_reasons.txt", sep="")
  }
  contig_ends <- read.table(data_file, header = FALSE, sep = ";")
  contig_end_data <- count(contig_ends$V1)
  names(contig_end_data) <- c("id", paste("Freq_", kmersize, sep = ""))
  
  for (kmersize in c(35, 55,75)) {
    if (reverse) {
      data_file <- paste("./", folder , "/results_for_",
                         kmersize, "/reverse_contig_end_reasons.txt", sep="")
    } else {
      data_file <- paste("./", folder , "/results_for_",
                         kmersize, "/contig_end_reasons.txt", sep="")
    }
    contig_ends <- read.table(data_file, header = FALSE, sep = ";")
    end_prop <- count(contig_ends$V1)
    names(end_prop) <- c("id", paste("Freq_", kmersize, sep = ""))
    contig_end_data <- merge(contig_end_data, end_prop, by.x="id", by.y="id")
  }
  
  contig_end_data <- data.frame(contig_end_data)
  contig_end_data <- melt(contig_end_data)
  abc <- ggplot(contig_end_data, aes(x = variable, y = value, fill = id )) +
    geom_bar( position = 'fill', colour="black", size=0.25, stat = 'identity' ) +
    coord_flip() + 
    labs(x = "K-mer Size") +
    scale_y_continuous(labels = percent_format()) +
    scale_x_discrete(labels = c("Freq_15" = "15","Freq_35" = "35",
                                "Freq_55" = "55","Freq_75" = "75")) +
    scale_fill_manual(labels = c("already visited, too small, x" = "Multiplicity\nDecreases         ",
                                 "already visited, too large, x" = "Multiplicity\nIncreases         ",
                                 "ambigiuous result, x" = "Ambigious\nBase              ",
                                 "too large, x" = "SI\nToo Large      ",
                                 "too many read ends, x" = "SI\nUnder-Populated   ",
                                 "too small, x" = "SI\nToo Small      ",
                                 "already visited, x" = "Already\nVisited             "),
                      values=c(
                        "already visited, too small, x" = "#0A369D",
                        "already visited, too large, x" = "#BBBE64",
                        "ambigiuous result, x" = "#731963",
                        "too large, x" = "#C3C3E6",
                        "too many read ends, x" = "#0B3C49",
                        "too small, x" = "#C1A5A9",
                        "already visited, x"="#CBD2D0"
                      )) +
    theme_tufte() + theme(axis.line = element_line(size = 0.25, colour = "black"),
                          axis.text.x = element_text(colour = "black"),
                          panel.grid.major=element_line(size = 0.1, colour = "gray"),
                          panel.grid.minor=element_line(size = 0.05, colour = "gray"),
                          legend.position="none", axis.title.x = element_blank())
    if (reverse) {
      plot_file <- file.path(paste("../../../WriteUp/plots/", folder,
                                   "/reverse_contig_end_reasons.png", sep=""))
    } else {
      plot_file <- file.path(paste("../../../WriteUp/plots/", folder,
                                   "/contig_end_reasons.png", sep=""))
    }
    ggsave(filename = plot_file, abc, width = 22, height = 4.5, units = "cm")
}



get_end_reasons_legends <- function() {

  contig_ends <- read.table("./legend_gen.txt", header = FALSE, sep = ";")
  contig_end_data <- count(contig_ends$V1)
  names(contig_end_data) <- c("id", paste("Freq_", 15, sep = ""))
  
  
  contig_end_data <- data.frame(contig_end_data)
  contig_end_data <- contig_end_data[1:4,]
  contig_end_data <- melt(contig_end_data)
  abc <- ggplot(contig_end_data, aes(x = variable, y = value, fill = id )) +
    geom_bar( position = 'fill', colour="black", size=0.25, stat = 'identity' ) +
    coord_flip() + 
    scale_fill_manual(labels = c("already visited, too small, x" = "Multiplicity\nDecreases            ",
                                 "already visited, too large, x" = "Multiplicity\nIncreases            ",
                                 "ambigiuous result, x" = "Ambigious\nBase                 ",
                                 "too large, x" = "Suffix Interval\nToo Large           ",
                                 "too many read ends, x" = "Suffix Interval\nUnder-Populated        ",
                                 "too small, x" = "Suffix Interval\nToo Small           ",
                                 "already visited, x" = "Already\nVisited        "),
                      values=c(
                        "already visited, too small, x" = "#0A369D",
                        "already visited, too large, x" = "#BBBE64",
                        "ambigiuous result, x" = "#731963",
                        "too large, x" = "#C3C3E6",
                        "too many read ends, x" = "#0B3C49",
                        "too small, x" = "#C1A5A9",
                        "already visited, x"="#CBD2D0"
                      )) +
    theme_tufte() + theme(axis.line = element_line(size = 0.25, colour = "black"),
                          axis.text.x = element_text(colour = "black"),
                          panel.grid.major=element_line(size = 0.1, colour = "gray"),
                          panel.grid.minor=element_line(size = 0.05, colour = "gray"),
                          legend.position="bottom", legend.title=element_blank(),
                          axis.title.x = element_blank()) + guides(guide_legend())
    plot_file <- file.path(paste("../../../WriteUp/plots/",
                                 "legend_for_contig_end_reasons_top.png", sep=""))
    ggsave(filename = plot_file, abc, width = 22, height = 5, units = "cm")
    
    
    contig_ends <- read.table("./legend_gen.txt", header = FALSE, sep = ";")
    contig_end_data <- count(contig_ends$V1)
    names(contig_end_data) <- c("id", paste("Freq_", 15, sep = ""))
    
    contig_end_data <- data.frame(contig_end_data)
    contig_end_data <- contig_end_data[5:7,]
    
    contig_end_data <- melt(contig_end_data)
    abc <- ggplot(contig_end_data, aes(x = variable, y = value, fill = id )) +
      geom_bar( position = 'fill', colour="black", size=0.25, stat = 'identity' ) +
      coord_flip() + 
      scale_fill_manual(labels=c("already visited, too small, x" = "Multiplicity\nDecreases           ",
                                "already visited, too large, x" = "Multiplicity\nIncreases           ",
                                 "ambigiuous result, x" = "Ambigious\nBase                  ",
                                 "too large, x" = "SI\nToo Large          ",
                                 "too many read ends, x" = "SI\nUnder-Populated        ",
                                 "too small, x" = "SI\nToo Small          ",
                                 "already visited, x" = "Already\nVisited          "),
                        values=c(
                          "already visited, too small, x" = "#0A369D",
                          "already visited, too large, x" = "#BBBE64",
                          "ambigiuous result, x" = "#731963",
                          "too large, x" = "#C3C3E6",
                          "too many read ends, x" = "#0B3C49",
                          "too small, x" = "#C1A5A9",
                          "already visited, x"="#CBD2D0"
                        )) +
      theme_tufte() + theme(axis.line = element_line(size = 0.25, colour = "black"),
                            axis.text.x = element_text(colour = "black"),
                            panel.grid.major=element_line(size = 0.1, colour = "gray"),
                            panel.grid.minor=element_line(size = 0.05, colour = "gray"),
                            legend.position="bottom", legend.title=element_blank(),
                            axis.title.x = element_blank()) + guides(guide_legend())
    plot_file <- file.path(paste("../../../WriteUp/plots/",
                                 "legend_for_contig_end_reasons_bottom.png", sep=""))
    
    ggsave(filename = plot_file, abc, width = 22, height = 5, units = "cm")
}



assembly_stats <- function(folder, kmersize) {
  # Take the sum of the contig lengths
  data_file <- paste("./", folder , "/results_for_", kmersize, "/contig_length.txt", sep="")
  contig_len <- read.table(data_file, header = FALSE, sep = ",", col.names=c("Length"))
  contig_len <- subset(contig_len, Length > kmersize*3)
  
  # full genome len:
  if (folder == "human30") { genome_len = 2237201 }
  else if (folder == "human60") { genome_len = 2237250 }
  else if (folder == "ERR431029") { genome_len = 6726781}
  else if (folder == "ERR430993") { genome_len = 7028270}
  
  print(paste("folder:", folder, " kmersize:", kmersize))
  print(paste("Genome Len:", genome_len, "Genome coverage:", sum(contig_len$Length)/genome_len))
  
  # N50
  sorted_len <- rev( sort(contig_len$Length) ) 
  cum_sum_bools <- cumsum(sorted_len) <= sum(contig_len$Length)/2 
  print(paste("N50:", sorted_len[ sum(cum_sum_bools) ], "L50:", sum(cum_sum_bools))) 
  
  # NG50
  sorted_len <- rev( sort(contig_len$Length) ) 
  cum_sum_bools <- cumsum(sorted_len) <= genome_len/2 
  print(paste("NG50:", sorted_len[ sum(cum_sum_bools) ], "LG50:", sum(cum_sum_bools))) 
  
  # Largest and average contig size:
  print(paste("Num Contigs:", nrow(contig_len), " Largest Contig:",  max(contig_len$Length), 
              " Avg Len:", mean(contig_len$Length) ))
  
  print(" ")
}




