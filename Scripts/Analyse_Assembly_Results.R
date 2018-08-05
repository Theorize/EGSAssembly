

source("/Users/isla/Documents/University/Masters/Thesis/Code/truncated-suffix-array/scripts/Analysis_Functions.R")

#####  BASE DISTRIBUTION 
## TO DO: same axes in every plot.
for (folder in c("human30", "human60", "ERR431029", "ERR430993")) {
    for (kmersize in c(15, 35, 55, 75)) {
    data_file <- paste("./", folder , "/results_for_", kmersize, "/base_dist.txt", sep="")
    plot_file <- file.path(paste("../../../WriteUp/plots/", folder, "/base_dist_", kmersize, ".png", sep=""))
    plot_data_dist(data_file, plot_file, FALSE)

    print(folder)
    print(kmersize)
    print('----------')
  }
}



###### CONTIG LENGTHS
create_contig_len_hist_human(5, "human30", 15)
create_contig_len_hist_human(50, "human30", 35)
create_contig_len_hist_human(200, "human30", 55)
create_contig_len_hist_human(200, "human30", 75)

create_contig_len_hist_human(5, "human60", 15)
create_contig_len_hist_human(50, "human60", 35)
create_contig_len_hist_human(200, "human60", 55)
create_contig_len_hist_human(200, "human60", 75)

create_contig_len_hist(1, "ERR431029", 15)
create_contig_len_hist(5, "ERR431029", 35)
create_contig_len_hist(10, "ERR431029", 55)
create_contig_len_hist(20, "ERR431029", 75)

create_contig_len_hist(1, "ERR430993", 15)
create_contig_len_hist(5, "ERR430993", 35)
create_contig_len_hist(10, "ERR430993", 55)
create_contig_len_hist(20, "ERR430993", 75)


plot_contig_end_reasons("human60", TRUE)
plot_contig_end_reasons("human30", TRUE)
plot_contig_end_reasons("human60", FALSE)
plot_contig_end_reasons("human30", FALSE)

plot_contig_end_reasons("ERR430993", FALSE)
plot_contig_end_reasons("ERR431029", FALSE)
get_end_reasons_legends()


for (folder in c("human30", "human60", "ERR431029", "ERR430993")) {  
  for (kmersize in c(15, 35, 55, 75)) {
    assembly_stats(folder, kmersize)
  }
}







###  MESSY CODE

############# For each k-mer size, plot a kmer spectrum.
setwd("/Users/isla/Documents/University/Masters/Thesis/Code/truncated-suffix-array/DataResults/ERR430993/kmer_stats")
for (kmersize in c(15, 35, 55, 75)) {
    data_file <- paste("./hist", kmersize, ".txt", sep="")
    if (file.exists(data_file)) {
        plot_file <- file.path(paste("./plots-hist", kmersize, ".png", sep=""))
        mydata <- read.table(data_file)
        mydata2<-mydata[0:400, ]
        mydata2 <- subset(mydata2, V2>0)
        
        point <- format_format(scientific = FALSE)
        
        abc <- qplot(V1,V2,data=mydata2, geom="line", main = paste("K-mer Size:", kmersize, sep = " "), xlab = "K-mer Frequency (Copy Number)", ylab = "Count of Distinct K-mers") + geom_point(size=0.75)
        #abc <- abc + ggtitle(paste("K-mer Length:", kmersize, " Expected K-mer Coverage:", expected_coverage, sep=" "))
        abc <- abc + scale_x_continuous(breaks = seq(0, 400, 50), minor_breaks = NULL)
        abc <- abc + scale_y_continuous(breaks = seq(0, 250000, 50000), labels = point, limits = c(0, 260000), minor_breaks = NULL)
        abc <- abc + theme_tufte()      # theme_tufte or theme_few      base_size   base_family
        abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major=element_line(size = 0.1, colour = "gray"), plot.title = element_text(hjust = 0.5))
        ggsave(filename = plot_file, abc)
        abc
        
        plot_file <- file.path(paste("./plots-cumulative", kmersize, ".png", sep=""))
        mydata2[,"cum_freq"] <- cumsum(mydata2$V2)
        abc <- qplot(V1,cum_freq,data=mydata2, geom="line", main = paste("K-mer Size:", kmersize, sep = " "), xlab = "K-mer Frequency (Copy Number)", ylab = "Count of Distinct K-mers") + geom_point(size=0.75)
        abc <- abc + scale_x_continuous(breaks = seq(0, 400, 50), minor_breaks = NULL)
        abc <- abc + scale_y_continuous(breaks = seq(0, 7000000, 1000000), labels = point, limits = c(0, 7250000), minor_breaks = NULL)
        abc <- abc + theme_tufte()      # theme_tufte or theme_few      base_size   base_family
        abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major=element_line(size = 0.1, colour = "gray"), plot.title = element_text(hjust = 0.5))
        ggsave(filename = plot_file, abc)
        abc
    }
}




setwd("/Users/isla/Documents/University/Masters/Thesis/Code/truncated-suffix-array/DataResults/ERR431029/kmer_stats")
for (kmersize in c(15, 35, 55, 75)) {
    data_file <- paste("./hist", kmersize, ".txt", sep="")
    if (file.exists(data_file)) {
        plot_file <- file.path(paste("./plots-hist", kmersize, ".png", sep=""))
        mydata <- read.table(data_file)
        mydata2<-mydata[0:400, ]
        mydata2 <- subset(mydata2, V2>0)
        
        point <- format_format(scientific = FALSE)
        
        abc <- qplot(V1,V2,data=mydata2, geom="line", main = paste("K-mer Size:", kmersize, sep = " "), xlab = "K-mer Frequency (Copy Number)", ylab = "Count of Distinct K-mers") + geom_point(size=0.75)
        abc <- abc + scale_x_continuous(breaks = seq(0, 400, 50), minor_breaks = NULL)
        abc <- abc + scale_y_continuous(breaks = seq(0, 350000, 50000), labels = point, limits = c(0, 350000), minor_breaks = NULL)
        abc <- abc + theme_tufte()      # theme_tufte or theme_few      base_size   base_family
        abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major=element_line(size = 0.1, colour = "gray"), plot.title = element_text(hjust = 0.5))
        ggsave(filename = plot_file, abc)
        abc
        
        
        plot_file <- file.path(paste("./plots-cumulative", kmersize, ".png", sep=""))
        mydata2[,"cum_freq"] <- cumsum(mydata2$V2)
        abc <- qplot(V1,cum_freq,data=mydata2, geom="line", main = paste("K-mer Size:", kmersize, sep = " "), xlab = "K-mer Frequency (Copy Number)", ylab = "Count of Distinct K-mers") + geom_point(size=0.75)
        abc <- abc + scale_x_continuous(breaks = seq(0, 400, 50), minor_breaks = NULL)
        abc <- abc + scale_y_continuous(breaks = seq(0, 7000000, 1000000), labels = point, limits = c(0, 7250000), minor_breaks = NULL)
        abc <- abc + theme_tufte()      # theme_tufte or theme_few      base_size   base_family
        abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major=element_line(size = 0.1, colour = "gray"), plot.title = element_text(hjust = 0.5))
        ggsave(filename = plot_file, abc, width = 11.5, height = 11.5, units = "cm")
        abc
    }
}






setwd("/Users/isla/Documents/University/Masters/Thesis/Code/truncated-suffix-array/DataResults/human30/kmer_stats")
# expected_coverage <- coverage * (read_length - kmersize + 1) / read_length

############# For each k-mer size, plot a kmer spectrum.
for (kmersize in c(15, 35, 55, 75)) {
    data_file <- paste("./hist", kmersize, ".txt", sep="")
    if (file.exists(data_file)) {
        plot_file <- file.path(paste("./plots-hist", kmersize, ".png", sep=""))
        mydata <- read.table(data_file)
        mydata2<-mydata[0:400, ]
        mydata2 <- subset(mydata2, V2>0)
        
        point <- format_format(scientific = FALSE)
        
        abc <- qplot(V1,V2,data=mydata2, geom="line", main = paste("K-mer Size:", kmersize, sep = " "), xlab = "K-mer Frequency (Copy Number)", ylab = "Count of Distinct K-mers") + geom_point(size=0.75)
        #abc <- abc + ggtitle(paste("K-mer Length:", kmersize, " Expected K-mer Coverage:", expected_coverage, sep=" "))
        abc <- abc + scale_x_continuous(breaks = seq(0, 400, 50), minor_breaks = NULL)
        abc <- abc + scale_y_continuous(breaks = seq(0, 800000, 100000), labels = point, limits = c(0, 810000), minor_breaks = NULL)
        abc <- abc + theme_tufte()      # theme_tufte or theme_few      base_size   base_family
        abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major=element_line(size = 0.1, colour = "gray"), plot.title = element_text(hjust = 0.5))
        ggsave(filename = plot_file, abc, width = 11.5, height = 11.5, units = "cm")
        abc
        
        plot_file <- file.path(paste("./plots-cumulative", kmersize, ".png", sep=""))
        mydata2[,"cum_freq"] <- cumsum(mydata2$V2)
        abc <- qplot(V1,cum_freq,data=mydata2, geom="line", main = paste("K-mer Size:", kmersize, sep = " "), xlab = "K-mer Frequency (Copy Number)", ylab = "Count of Distinct K-mers") + geom_point(size=0.75)
        abc <- abc + scale_x_continuous(breaks = seq(0, 400, 50), minor_breaks = NULL)
        abc <- abc + scale_y_continuous(breaks = seq(0, 2000000, 250000), labels = point, limits = c(0, 2000000), minor_breaks = NULL)
        abc <- abc + theme_tufte()      # theme_tufte or theme_few      base_size   base_family
        abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major=element_line(size = 0.1, colour = "gray"), plot.title = element_text(hjust = 0.5))
        ggsave(filename = plot_file, abc, width = 11.5, height = 11.5, units = "cm")
        abc
    }
}










setwd("/Users/isla/Documents/University/Masters/Thesis/Code/truncated-suffix-array/DataResults/human60/kmer_stats")
# expected_coverage <- coverage * (read_length - kmersize + 1) / read_length

############# For each k-mer size, plot a kmer spectrum.
for (kmersize in c(15, 35, 55, 75)) {
    data_file <- paste("./hist", kmersize, ".txt", sep="")
    if (file.exists(data_file)) {
        plot_file <- file.path(paste("./plots-hist", kmersize, ".png", sep=""))
        mydata <- read.table(data_file)
        mydata2<-mydata[0:400, ]
        mydata2 <- subset(mydata2, V2>0)
        
        point <- format_format(scientific = FALSE)
        
        abc <- qplot(V1,V2,data=mydata2, geom="line", main = paste("K-mer Size:", kmersize, sep = " "), xlab = "K-mer Frequency (Copy Number)", ylab = "Count of Distinct K-mers") + geom_point(size=0.75)
        #abc <- abc + ggtitle(paste("K-mer Length:", kmersize, " Expected K-mer Coverage:", expected_coverage, sep=" "))
        abc <- abc + scale_x_continuous(breaks = seq(0, 400, 50), minor_breaks = NULL)
        abc <- abc + scale_y_continuous(breaks = seq(0, 1750000, 250000), labels = point, limits = c(0, 1700000), minor_breaks = NULL)
        abc <- abc + theme_tufte()      # theme_tufte or theme_few      base_size   base_family
        abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major=element_line(size = 0.1, colour = "gray"), plot.title = element_text(hjust = 0.5))
        ggsave(filename = plot_file, abc)
        abc
        
        plot_file <- file.path(paste("./plots-cumulative", kmersize, ".png", sep=""))
        mydata2[,"cum_freq"] <- cumsum(mydata2$V2)
        abc <- qplot(V1,cum_freq,data=mydata2, geom="line", main = paste("K-mer Size:", kmersize, sep = " "), xlab = "K-mer Frequency (Copy Number)", ylab = "Count of Distinct K-mers") + geom_point(size=0.75)
        abc <- abc + scale_x_continuous(breaks = seq(0, 400, 50), minor_breaks = NULL)
        abc <- abc + scale_y_continuous(breaks = seq(0, 2000000, 250000), labels = point, limits = c(0, 2000000), minor_breaks = NULL)
        abc <- abc + theme_tufte()      # theme_tufte or theme_few      base_size   base_family
        abc <- abc + theme(axis.line = element_line(size = 0.25, colour = "black"), panel.grid.major=element_line(size = 0.1, colour = "gray"), plot.title = element_text(hjust = 0.5))
        ggsave(filename = plot_file, abc)
        abc
    }
}






