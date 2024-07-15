if (!require("ggpubr")) {
  install.packages("ggpubr", dependencies = TRUE)
  library(ggpubr)}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)}
if (!require("pacman")) {
  install.packages("pacman", dependencies = TRUE)
  library(pacman)}
if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(pacman)}

####File preparation####

#Load your table of GT, names and positions. Thhe input files are achieved thorugh a 012 transformation from a vcf tfile using vcftools --012
table <-read.table(file = "chr_biallele_005maf_5DP_Raw.012", 
                   sep = "\t", header=FALSE)
names <- read.table(file = "chr_biallele_005maf_5DP_Raw.012.indv", 
                   sep = "\t", header=FALSE)
position <- read.table(file = "chr_biallele_005maf_5DP_Raw.012.pos", 
                      sep = "\t", header=FALSE) 

table <- as.data.frame(table)
names <- as.data.frame(names)
position <- as.data.frame(position)

#Add headers and remove unnecessary columns
positionhead <- c("CHR", "POS")
colnames(position) <- positionhead
table <- select(table, c(-1))

#Transpose the matrices and get the input file for the Detect_Inversions function
genotips <- t(table)
genotips <- as.data.frame(genotips)
noms <- t(names)
noms <- as.data.frame(names)
position <- as.data.frame(position)
colnames(genotips) <- names
TS <- cbind(position,genotips)
TS[TS == -1] <- NA


####Detect_Inversions Function####

Detect_Inversions <- function(df, window, step) {
  chromosome <- df[, 1]
  data_cols <- df[, -1]  # All columns except the first (chromosome)
  total <- nrow(df)
  spots <- seq(from = 1, to = (total - window + 1), by = step)
  
  # Initialize result data frame with columns for each data column's mean and chromosome info
  result <- data.frame(
    chromosome = character(length(spots)),
    stringsAsFactors = FALSE
  )
  
  for (col in names(data_cols)) {
    result[[paste0(col, "_mean")]] <- numeric(length(spots))
  }
  
  for (i in 1:length(spots)) {
    window_chromosome <- chromosome[spots[i]:(spots[i] + window - 1)]
    if (length(unique(window_chromosome)) == 1) {  # All chromosomes in the window are the same
      result$chromosome[i] <- unique(window_chromosome)
      for (col in names(data_cols)) {
        window_data <- data_cols[[col]][spots[i]:(spots[i] + window - 1)]
        if (all(is.na(window_data))) {  # If all values in the window are NA
          result[[paste0(col, "_mean")]][i] <- NA
        } else {
          result[[paste0(col, "_mean")]][i] <- mean(window_data, na.rm = TRUE)  # Calculate mean excluding NAs
        }
      }
    } else {
      result$chromosome[i] <- window_chromosome[1]  # Keep the first chromosome in the window
      for (col in names(data_cols)) {
        window_data <- data_cols[[col]][spots[i]:(spots[i] + window - 1)]
        if (all(is.na(window_data))) {  # If all values in the window are NA
          result[[paste0(col, "_mean")]][i] <- NA
        } else {
          result[[paste0(col, "_mean")]][i] <- mean(window_data, na.rm = TRUE)  # Calculate mean excluding NAs
        }
      }
    }
  }
  
  return(result)
}

####Usage####

newTS <- Detect_Inversions(TS, window=10000, step = 2500)
newTS <- as.data.frame(newTS)

library(reshape2)
newTS$POS_mean <- as.character(newTS$POS_mean)
newTS<-melt(newTS)
newTS$POS_mean <- as.numeric(newTS$POS_mean)


pdf(file="inversion_detection.pdf", width=10, height=15)
ggplot(data=newTS, aes(x=POS_mean, y=value, fill=variable, color=pop)) +
  geom_line(alpha=0.4)+
  theme_classic()+
  facet_grid(chromosome~., scales="free_y")+
  scale_x_continuous(labels = scales::comma)+                                                         
  scale_y_continuous(name="Homozygosity")+
  #scale_color_manual(#set your data here) +
  theme(legend.position = "bottom",
        axis.text=element_text(size=14, face="bold"),
        axis.title=element_text(size=20,face="bold"), axis.text.x =element_text(angle=0, hjust=1)) 
dev.off()
