
#Install the packages below if you dont already have them. Some may be through bioconductor
#Library Packages
library(DiffBind)
library(profileplyr)
library(dplyr)
library(ggplot2)
library(reshape)
library(stringr)
library(rtracklayer)

if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
library(tidyr)


#Read file that specifies file locations and information for DiffBind
#To keep color schemes consistent later, highly recommend putting your mutant / comparing construct 2nd in the data sheet, and wt as the 3rd
#This should be the only file you need to make changes, change sample names in the file and paths to the appropriate files
#then changing the file path below to where the "files.csv" is saved to
data_files <- read.csv("/Users/ayh/Downloads/Pipeline_PHBvWtD/outputs/files.csv", sep = ",", header = TRUE)

# Read the files.csv
file_paths <- read.csv("./outputs/files.csv", sep = ",", header = TRUE)$file_column_name

# Replace 'file_column_name' with the name of the column in your CSV that contains the file paths.

# Check if each file exists
file_check <- sapply(file_paths, file.exists)

# Print the results
if (all(file_check)) {
  cat("All referenced files exist and are accessible.\n")
} else {
  cat("Some files are missing or inaccessible. Check the following:\n")
  missing_files <- file_paths[!file_check]
  print(missing_files)
}
# Perform DiffBind Analysis
# Common error seen with dba.count: "Read block operation failed with error 1 after 0 of 4 bytes"
# Online says this error occurs when the bam file is corrupted during download and should be redownloaded. 
# This code is set up for an analysis of 2 experimental genotypes, with 3 replicates each, and the wildtype control with 2 replicates. If your analysis is different from that, you will need to adjust some of the code 
Comparison <- dba(sampleSheet= data_files)
Comparison <- dba.count(Comparison, summits = 150, filter = 0, minOverlap = 2)
Comparison <- dba.normalize(Comparison)
Comparison <- dba.contrast(Comparison, minMembers = 2)
Comparison <- dba.analyze(Comparison)

# Generate the report of the DiffBind analysis, convert to data frame, arrange in same order
ChIP.DB.1 <- dba.report(Comparison, th = 1, bNormalized = TRUE, bCalled = TRUE, contrast = 1)
ChIP.DB.2 <- dba.report(Comparison, th = 1, bNormalized = TRUE, bCalled = TRUE, contrast = 2)
ChIP.DB.3 <- dba.report(Comparison, th = 1, bNormalized = TRUE, bCalled = TRUE, contrast = 3)

ChIP.DB.1 <- as.data.frame(ChIP.DB.1)
ChIP.DB.1$row_name <- as.numeric(row.names(ChIP.DB.1))
ChIP.DB.1 <- ChIP.DB.1[order(ChIP.DB.1$row_name), ]

ChIP.DB.2 <- as.data.frame(ChIP.DB.2)
ChIP.DB.2$row_name <- as.numeric(row.names(ChIP.DB.2))
ChIP.DB.2 <- ChIP.DB.2[order(ChIP.DB.2$row_name), ]

ChIP.DB.3 <- as.data.frame(ChIP.DB.3)
ChIP.DB.3$row_name <- as.numeric(row.names(ChIP.DB.3))
ChIP.DB.3 <- ChIP.DB.3[order(ChIP.DB.3$row_name), ]

# Combine information from the data sets
db <- data.frame(ChIP.DB.1[1:5], ChIP.DB.1[7:8], ChIP.DB.2[8], ChIP.DB.1[11], ChIP.DB.2[11], ChIP.DB.3[11], 
                 ChIP.DB.1[12:13], ChIP.DB.2[13], ChIP.DB.1[14])
db1 <- db  # Keep the original data for reference

# Identify sites we call as true binding, with at least 2 replicates in one genotype having the same called peak
db <- filter(db, db[[12]] >= 2 | db[[13]] >= 2 | db[[14]] >= 2)

# Initialize remaining peaks
remaining_peaks <- db


# Step 9: Process Peaks by Categories
assign_relationship <- function(data, condition, relationship_label) {
  subset <- filter(data, !!rlang::parse_expr(condition))
  if (nrow(subset) > 0) {
    subset$relationship <- relationship_label
    data <- anti_join(data, subset, by = "row_name")
  }
  list(data = data, subset = subset)
}

# Step 1: Prioritize Unique Binding Sites
result <- assign_relationship(remaining_peaks,
                              "Called1 >= 2 & Called2 == 0 & Called2.1 == 0 & FDR <= 0.05 & FDR.1 <= 0.05 & Conc_PHB > Conc_WtD & Conc_PHB > Conc_Wt",
                              "sample1 unique"
)
remaining_peaks <- result$data
sample1_unique <- result$subset

result <- assign_relationship(remaining_peaks,
                              "Called1 == 0 & Called2 >= 2 & Called2.1 == 0 & FDR <= 0.05 & FDR.2 <= 0.05 & Conc_WtD > Conc_PHB & Conc_WtD > Conc_Wt",
                              "sample2 unique"
)
remaining_peaks <- result$data
sample2_unique <- result$subset

result <- assign_relationship(remaining_peaks,
                              "Called1 == 0 & Called2 == 0 & Called2.1 >= 2 & FDR.1 <= 0.05 & FDR.2 <= 0.05 & Conc_Wt > Conc_PHB & Conc_Wt > Conc_WtD",
                              "sample3 unique"
)
remaining_peaks <- result$data
sample3_unique <- result$subset

# Step 2: Identify True Mutuals
result <- assign_relationship(remaining_peaks,
                              "Called1 >= 1 & Called2 >= 1 & Called2.1 >= 1 & FDR >= 0.1 & FDR.1 >= 0.1 & FDR.2 >= 0.1",
                              "true_mutual"
)
remaining_peaks <- result$data
true_mutual <- result$subset

# Step 3: Refine Mutuals by Specific Categories
categories <- list(
  list(condition = "Called1 >= 1 & Called2 >= 1 & Called2.1 == 0 & FDR <= 0.05 & FDR.1 <= 0.05 & Conc_PHB > Conc_WtD",
       label = "mutual of samples1&2; sample1 high"),
  list(condition = "Called1 >= 1 & Called2 >= 1 & Called2.1 == 0 & FDR <= 0.05 & FDR.1 <= 0.05 & Conc_WtD > Conc_PHB",
       label = "mutual of samples1&2; sample2 high"),
  list(condition = "Called1 >= 1 & Called2 == 0 & Called2.1 >= 1 & FDR <= 0.05 & FDR.2 <= 0.05 & Conc_PHB > Conc_Wt",
       label = "mutual of samples1&3; sample1 high"),
  list(condition = "Called1 >= 1 & Called2 == 0 & Called2.1 >= 1 & FDR <= 0.05 & FDR.2 <= 0.05 & Conc_Wt > Conc_PHB",
       label = "mutual of samples1&3; sample3 high"),
  list(condition = "Called1 == 0 & Called2 >= 1 & Called2.1 >= 1 & FDR.1 <= 0.05 & FDR.2 <= 0.05 & Conc_WtD > Conc_Wt",
       label = "mutual of samples2&3; sample2 high"),
  list(condition = "Called1 == 0 & Called2 >= 1 & Called2.1 >= 1 & FDR.1 <= 0.05 & FDR.2 <= 0.05 & Conc_Wt > Conc_WtD",
       label = "mutual of samples2&3; sample3 high"),
  list(condition = "Called1 >= 1 & Called2 >= 1 & Called2.1 == 0 & FDR >= 0.05 & FDR.1 <= 0.05 & FDR.2 <= 0.05",
       label = "mutual of samples1&2"),
  list(condition = "Called1 >= 1 & Called2 == 0 & Called2.1 >= 1 & FDR <= 0.05 & FDR.1 >= 0.05 & FDR.2 <= 0.05",
       label = "mutual of samples1&3"),
  list(condition = "Called1 == 0 & Called2 >= 1 & Called2.1 >= 1 & FDR <= 0.05 & FDR.1 <= 0.05 & FDR.2 >= 0.05",
       label = "mutual of samples2&3")
)

for (cat in categories) {
  result <- assign_relationship(remaining_peaks, cat$condition, cat$label)
  remaining_peaks <- result$data
  assign(cat$label, result$subset)
}

# Classify Remaining Peaks as Removed
removed <- remaining_peaks
removed$relationship <- "removed"

# Combine all relationships
combined_df <- bind_rows(true_mutual, sample1_unique, sample2_unique, sample3_unique, 
                         `mutual of samples1&2; sample1 high`, `mutual of samples1&2; sample2 high`, 
                         `mutual of samples1&3; sample1 high`, `mutual of samples1&3; sample3 high`,
                         `mutual of samples2&3; sample2 high`, `mutual of samples2&3; sample3 high`,
                         `mutual of samples1&2`, `mutual of samples1&3`, `mutual of samples2&3`, removed)

# Write the output to CSV
combined_df <- combined_df[order(combined_df$row_name), ]
write.csv(combined_df, "diffbind_output_refined.csv", row.names = FALSE)



# Generate a summary of relationships
summary_df <- combined_df %>%
  group_by(relationship) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = relationship, values_from = count, values_fill = 0)
write.csv(summary_df, "relationship_summary_refined.csv", row.names = FALSE)


###########################################################
#Begin to Generate Heatmaps
#This code assumes that "sample3" is wildtype
#There are 2 codes that can alter the color (scale_fill_gradiant & scale_fill_distiller)
#Distiller will make the scale itself for the most part except the NA color which you should assign as the darkest color
#Depending on your level of signal, you may need to adjust the "limits = c(0,150) to better highlight your own data. 
#Gradiant is for if you want to make your own color color

#Generate heatmaps of the mutual sites of 1&2, mutual affinity, removing wt of sample 3
#Arrange data by FDR and Get the signal information
mutual_of_1_2 <- mutual_of_1_2[order(mutual_of_1_2$FDR), ]
sites <- mutual_of_1_2[1:3]
sites$row_name <- mutual_of_1_2$row_name
coordinates <- GRanges(sites[1:3])
coordinates$row_name <- row.names(mutual_of_1_2)
profile <- dba.plotProfile(Comparison, merge=c(DBA_CONDITION, DBA_REPLICATE), sites = coordinates, maxSites = 15000, distanceAround = 1000, labels = DBA_CONDITION)
heat <- as.data.frame(convertToEnrichedHeatmapMat(profile))
heat$row_name <- as.double(profile@rowRanges$names)
heat <- left_join(sites, heat, by = "row_name")

#Isolate the signal information for each sample
sample1_heat <- heat[5:104]
sample2_heat <- heat[105:204]
sample3_heat <- heat[205:304]

#convert heat values to proper formats
sample1_heat <- melt(sample1_heat)
sample1_heat <- data.frame(sample1_heat, gene = 1:nrow(heat))

sample2_heat <- melt(sample2_heat)
sample2_heat <- data.frame(sample2_heat, gene = 1:nrow(heat))

sample3_heat <- melt(sample3_heat)
sample3_heat <- data.frame(sample3_heat, gene = 1:nrow(heat))

#Generate Heatmap of Sample 1
g1 <- ggplot(sample1_heat, aes(x = variable, y=gene, fill = value))
g1 <- g1 + geom_raster(interpolate = TRUE)
#g1 <- g1 + scale_fill_gradient(low = "#FFF5F0" , high = "#FF0000", na.value = "#FF0000", limits = c(0,150))
g1 <- g1 + scale_fill_distiller(palette = "Reds", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#A0101D")
g1 <- g1 + theme_classic()
g1 <- g1 + scale_y_reverse()
g1 <- g1 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g1
ggsave(plot = g1, filename = "mutual of samples 1&2, sample 1.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 2
g2 <- ggplot(sample2_heat, aes(x = variable, y=gene, fill = value))
g2 <- g2 + geom_raster(interpolate = TRUE)
#Use this code below if you want to use a custom color set instead of the distiller palettes.
#g2 <- g2 + scale_fill_gradient(low = "#FFF5F0" , high = "#CC0066", na.value = "#CC0066", limits = c(0,150))
g2 <- g2 + scale_fill_distiller(palette = "Blues", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#092F6B")
g2 <- g2 + theme_classic()
g2 <- g2 + scale_y_reverse()
g2 <- g2 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g2
ggsave(plot = g2, filename = "mutual of samples 1&2, sample 2.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 3
g3 <- ggplot(sample3_heat, aes(x = variable, y=gene, fill = value))
g3 <- g3 + geom_raster(interpolate = TRUE)
#g3 <- g3 + scale_fill_gradient(low = "#FFF5F0" , high = "#000000", na.value = "#000000", limits = c(0,150))
g3 <- g3 + scale_fill_distiller(palette = "Greys", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#000000")
g3 <- g3 + theme_classic()
g3 <- g3 + scale_y_reverse()
g3 <- g3 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g3
ggsave(plot = g3, filename = "mutual of samples 1&2, sample 3.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Averages of Heat Signals for each sample
Avg_sample_heat <- data.frame("sample1" = colMeans(heat[5:104]), "sample2" = colMeans(heat[105:204]), "sample3" = colMeans(heat[205:304]))
write.csv(Avg_sample_heat, "mutual of samples 1&2 averaged heat values.csv", row.names = FALSE)


#Generate heatmaps of the mutual sites of 1&2, sample 1 higher affinity, removing wt of sample 3
#Arrange data by FDR and Get the signal information
mutual_of_1_2_1high <- mutual_of_1_2_1high[order(mutual_of_1_2_1high$FDR), ]
sites <- mutual_of_1_2_1high[1:3]
sites$row_name <- mutual_of_1_2_1high$row_name
coordinates <- GRanges(sites[1:3])
coordinates$row_name <- row.names(mutual_of_1_2_1high)
profile <- dba.plotProfile(Comparison, merge=c(DBA_CONDITION, DBA_REPLICATE), sites = coordinates, maxSites = 15000, distanceAround = 1000, labels = DBA_CONDITION)
heat <- as.data.frame(convertToEnrichedHeatmapMat(profile))
heat$row_name <- as.double(profile@rowRanges$names)
heat <- left_join(sites, heat, by = "row_name")

#Isolate the signal information for each sample
sample1_heat <- heat[5:104]
sample2_heat <- heat[105:204]
sample3_heat <- heat[205:304]

#convert heat values to proper formats
sample1_heat <- melt(sample1_heat)
sample1_heat <- data.frame(sample1_heat, gene = 1:nrow(heat))

sample2_heat <- melt(sample2_heat)
sample2_heat <- data.frame(sample2_heat, gene = 1:nrow(heat))

sample3_heat <- melt(sample3_heat)
sample3_heat <- data.frame(sample3_heat, gene = 1:nrow(heat))

#Generate Heatmap of Sample 1
g1 <- ggplot(sample1_heat, aes(x = variable, y=gene, fill = value))
g1 <- g1 + geom_raster(interpolate = TRUE)
#g1 <- g1 + scale_fill_gradient(low = "#FFF5F0" , high = "#FF0000", na.value = "#FF0000", limits = c(0,150))
g1 <- g1 + scale_fill_distiller(palette = "Reds", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#A0101D")
g1 <- g1 + theme_classic()
g1 <- g1 + scale_y_reverse()
g1 <- g1 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g1
ggsave(plot = g1, filename = "mutual of samples 1&2, sample 1 high, sample 1.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 2
g2 <- ggplot(sample2_heat, aes(x = variable, y=gene, fill = value))
g2 <- g2 + geom_raster(interpolate = TRUE)
#g2 <- g2 + scale_fill_gradient(low = "#FFF5F0" , high = "#CC0066", na.value = "#CC0066", limits = c(0,150))
g2 <- g2 + scale_fill_distiller(palette = "Blues", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#092F6B")
g2 <- g2 + theme_classic()
g2 <- g2 + scale_y_reverse()
g2 <- g2 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g2
ggsave(plot = g2, filename = "mutual of samples 1&2, sample 1 high, sample 2.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 3
g3 <- ggplot(sample3_heat, aes(x = variable, y=gene, fill = value))
g3 <- g3 + geom_raster(interpolate = TRUE)
#g3 <- g3 + scale_fill_gradient(low = "#FFF5F0" , high = "#000000", na.value = "#000000", limits = c(0,150))
g3 <- g3 + scale_fill_distiller(palette = "Greys", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#000000")
g3 <- g3 + theme_classic()
g3 <- g3 + scale_y_reverse()
g3 <- g3 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g3
ggsave(plot = g3, filename = "mutual of samples 1&2, sample 1 high, sample 3.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Averages of Heat Signals for each sample
Avg_sample_heat <- data.frame("sample1" = colMeans(heat[5:104]), "sample2" = colMeans(heat[105:204]), "sample3" = colMeans(heat[205:304]))
write.csv(Avg_sample_heat, "mutual of samples 1&2, sample 1 high averaged heat values.csv", row.names = FALSE)


#Generate heatmaps of the mutual sites of 1&2, sample 2 higher affinity, removing wt of sample 3
#Arrange data by FDR and Get the signal information
mutual_of_1_2_2high <- mutual_of_1_2_2high[order(mutual_of_1_2_2high$FDR), ]
sites <- mutual_of_1_2_2high[1:3]
sites$row_name <- mutual_of_1_2_2high$row_name
coordinates <- GRanges(sites[1:3])
coordinates$row_name <- row.names(mutual_of_1_2_2high)
profile <- dba.plotProfile(Comparison, merge=c(DBA_CONDITION, DBA_REPLICATE), sites = coordinates, maxSites = 15000, distanceAround = 1000, labels = DBA_CONDITION)
heat <- as.data.frame(convertToEnrichedHeatmapMat(profile))
heat$row_name <- as.double(profile@rowRanges$names)
heat <- left_join(sites, heat, by = "row_name")

#Isolate the signal information for each sample
sample1_heat <- heat[5:104]
sample2_heat <- heat[105:204]
sample3_heat <- heat[205:304]

#convert heat values to proper formats
sample1_heat <- melt(sample1_heat)
sample1_heat <- data.frame(sample1_heat, gene = 1:nrow(heat))

sample2_heat <- melt(sample2_heat)
sample2_heat <- data.frame(sample2_heat, gene = 1:nrow(heat))

sample3_heat <- melt(sample3_heat)
sample3_heat <- data.frame(sample3_heat, gene = 1:nrow(heat))

#Generate Heatmap of Sample 1
g1 <- ggplot(sample1_heat, aes(x = variable, y=gene, fill = value))
g1 <- g1 + geom_raster(interpolate = TRUE)
#g1 <- g1 + scale_fill_gradient(low = "#FFF5F0" , high = "#FF0000", na.value = "#FF0000", limits = c(0,150))
g1 <- g1 + scale_fill_distiller(palette = "Reds", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#A0101D")
g1 <- g1 + theme_classic()
g1 <- g1 + scale_y_reverse()
g1 <- g1 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g1
ggsave(plot = g1, filename = "mutual of samples 1&2, sample 2 high, sample 1.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 2
g2 <- ggplot(sample2_heat, aes(x = variable, y=gene, fill = value))
g2 <- g2 + geom_raster(interpolate = TRUE)
#g2 <- g2 + scale_fill_gradient(low = "#FFF5F0" , high = "#CC0066", na.value = "#CC0066", limits = c(0,150))
g2 <- g2 + scale_fill_distiller(palette = "Blues", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#092F6B")
g2 <- g2 + theme_classic()
g2 <- g2 + scale_y_reverse()
g2 <- g2 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g2
ggsave(plot = g2, filename = "mutual of samples 1&2, sample 2 high, sample 2.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 3
g3 <- ggplot(sample3_heat, aes(x = variable, y=gene, fill = value))
g3 <- g3 + geom_raster(interpolate = TRUE)
#g3 <- g3 + scale_fill_gradient(low = "#FFF5F0" , high = "#000000", na.value = "#000000", limits = c(0,150))
g3 <- g3 + scale_fill_distiller(palette = "Greys", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#000000")
g3 <- g3 + theme_classic()
g3 <- g3 + scale_y_reverse()
g3 <- g3 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g3
ggsave(plot = g3, filename = "mutual of samples 1&2, sample 2 high, sample 3.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Averages of Heat Signals for each sample
Avg_sample_heat <- data.frame("sample1" = colMeans(heat[5:104]), "sample2" = colMeans(heat[105:204]), "sample3" = colMeans(heat[205:304]))
write.csv(Avg_sample_heat, "mutual of samples 1&2, sample 2 high averaged heat values.csv", row.names = FALSE)


#Generate heatmaps of the sample 1 unique
#Arrange data by FDR and Get the signal information
sample1_unique <- sample1_unique[order(sample1_unique$FDR), ]
sites <- sample1_unique[1:3]
sites$row_name <- sample1_unique$row_name
coordinates <- GRanges(sites[1:3])
coordinates$row_name <- row.names(sample1_unique)
profile <- dba.plotProfile(Comparison, merge=c(DBA_CONDITION, DBA_REPLICATE), sites = coordinates, maxSites = 15000, distanceAround = 1000, labels = DBA_CONDITION)
heat <- as.data.frame(convertToEnrichedHeatmapMat(profile))
heat$row_name <- as.double(profile@rowRanges$names)
heat <- left_join(sites, heat, by = "row_name")

#Isolate the signal information for each sample
sample1_heat <- heat[5:104]
sample2_heat <- heat[105:204]
sample3_heat <- heat[205:304]

#convert heat values to proper formats
sample1_heat <- melt(sample1_heat)
sample1_heat <- data.frame(sample1_heat, gene = 1:nrow(heat))

sample2_heat <- melt(sample2_heat)
sample2_heat <- data.frame(sample2_heat, gene = 1:nrow(heat))

sample3_heat <- melt(sample3_heat)
sample3_heat <- data.frame(sample3_heat, gene = 1:nrow(heat))

#Generate Heatmap of Sample 1
g1 <- ggplot(sample1_heat, aes(x = variable, y=gene, fill = value))
g1 <- g1 + geom_raster(interpolate = TRUE)
#g1 <- g1 + scale_fill_gradient(low = "#FFF5F0" , high = "#FF0000", na.value = "#FF0000", limits = c(0,150))
g1 <- g1 + scale_fill_distiller(palette = "Reds", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#A0101D")
g1 <- g1 + theme_classic()
g1 <- g1 + scale_y_reverse()
g1 <- g1 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g1
ggsave(plot = g1, filename = "sample 1 unique, sample 1.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 2
g2 <- ggplot(sample2_heat, aes(x = variable, y=gene, fill = value))
g2 <- g2 + geom_raster(interpolate = TRUE)
#g2 <- g2 + scale_fill_gradient(low = "#FFF5F0" , high = "#CC0066", na.value = "#CC0066", limits = c(0,150))
g2 <- g2 + scale_fill_distiller(palette = "Blues", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#092F6B")
g2 <- g2 + theme_classic()
g2 <- g2 + scale_y_reverse()
g2 <- g2 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g2
ggsave(plot = g2, filename = "sample 1 unique, sample 2.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 3
g3 <- ggplot(sample3_heat, aes(x = variable, y=gene, fill = value))
g3 <- g3 + geom_raster(interpolate = TRUE)
#g3 <- g3 + scale_fill_gradient(low = "#FFF5F0" , high = "#000000", na.value = "#000000", limits = c(0,150))
g3 <- g3 + scale_fill_distiller(palette = "Greys", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#000000")
g3 <- g3 + theme_classic()
g3 <- g3 + scale_y_reverse()
g3 <- g3 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g3
ggsave(plot = g3, filename = "sample 1 unique, sample 3.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Averages of Heat Signals for each sample
Avg_sample_heat <- data.frame("sample1" = colMeans(heat[5:104]), "sample2" = colMeans(heat[105:204]), "sample3" = colMeans(heat[205:304]))
write.csv(Avg_sample_heat, "sample 1 unique averaged heat values.csv", row.names = FALSE)


#Generate heatmaps of the sample 2 unique
#Arrange data by FDR and Get the signal information
sample2_unique <- sample2_unique[order(sample2_unique$FDR), ]
sites <- sample2_unique[1:3]
sites$row_name <- sample2_unique$row_name
coordinates <- GRanges(sites[1:3])
coordinates$row_name <- row.names(sample2_unique)
profile <- dba.plotProfile(Comparison, merge=c(DBA_CONDITION, DBA_REPLICATE), sites = coordinates, maxSites = 15000, distanceAround = 1000, labels = DBA_CONDITION)
heat <- as.data.frame(convertToEnrichedHeatmapMat(profile))
heat$row_name <- as.double(profile@rowRanges$names)
heat <- left_join(sites, heat, by = "row_name")

#Isolate the signal information for each sample
sample1_heat <- heat[5:104]
sample2_heat <- heat[105:204]
sample3_heat <- heat[205:304]

#convert heat values to proper formats
sample1_heat <- melt(sample1_heat)
sample1_heat <- data.frame(sample1_heat, gene = 1:nrow(heat))

sample2_heat <- melt(sample2_heat)
sample2_heat <- data.frame(sample2_heat, gene = 1:nrow(heat))

sample3_heat <- melt(sample3_heat)
sample3_heat <- data.frame(sample3_heat, gene = 1:nrow(heat))

#Generate Heatmap of Sample 1
g1 <- ggplot(sample1_heat, aes(x = variable, y=gene, fill = value))
g1 <- g1 + geom_raster(interpolate = TRUE)
#g1 <- g1 + scale_fill_gradient(low = "#FFF5F0" , high = "#FF0000", na.value = "#FF0000", limits = c(0,150))
g1 <- g1 + scale_fill_distiller(palette = "Reds", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#A0101D")
g1 <- g1 + theme_classic()
g1 <- g1 + scale_y_reverse()
g1 <- g1 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g1
ggsave(plot = g1, filename = "sample 2 unique, sample 1.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 2
g2 <- ggplot(sample2_heat, aes(x = variable, y=gene, fill = value))
g2 <- g2 + geom_raster(interpolate = TRUE)
g2 <- g2 + scale_fill_gradient(low = "#FFF5F0" , high = "#CC0066", na.value = "#CC0066", limits = c(0,150))
#g2 <- g2 + scale_fill_distiller(palette = "Blues", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#092F6B")
g2 <- g2 + theme_classic()
g2 <- g2 + scale_y_reverse()
g2 <- g2 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g2
ggsave(plot = g2, filename = "sample 2 unique, sample 2.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Heatmap of Sample 3
g3 <- ggplot(sample3_heat, aes(x = variable, y=gene, fill = value))
g3 <- g3 + geom_raster(interpolate = TRUE)
#g3 <- g3 + scale_fill_gradient(low = "#FFF5F0" , high = "#000000", na.value = "#000000", limits = c(0,150))
g3 <- g3 + scale_fill_distiller(palette = "Greys", direction = 1, aesthetics = "fill", limits = c(0,150),na.value = "#000000")
g3 <- g3 + theme_classic()
g3 <- g3 + scale_y_reverse()
g3 <- g3 + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank(),
                 axis.text.y=element_blank(),
                 axis.ticks.y=element_blank(),
                 axis.title =element_blank(),
                 axis.line =element_blank(),
                 legend.position = "none",
                 panel.background = element_rect(fill='transparent'),
                 plot.background = element_rect(fill='transparent', color=NA),
)
g3
ggsave(plot = g3, filename = "sample 2 unique, sample 3.pdf", height = 1+(nrow(heat)/1000), width = 1, unit = "in")

#Generate Averages of Heat Signals for each sample
Avg_sample_heat <- data.frame("sample1" = colMeans(heat[5:104]), "sample2" = colMeans(heat[105:204]), "sample3" = colMeans(heat[205:304]))
write.csv(Avg_sample_heat, "sample 2 unique averaged heat values.csv", row.names = FALSE)


###################################################
#Begin ChIP Peak Annotation
# Define the path to your GFF3 file
# Easiest to just put it in the folder where the R project and this code is
gff3_file <- "Arabidopsis_thaliana.TAIR10.56.gff3"

# Read the GFF3 file
gff3_data <- import(gff3_file)
gff3_data <- as.data.frame(gff3_data)

# Filter the GFF3 data for "gene" and "ncRNA_gene" and extract relevant columns
filtered_gff3_data <- gff3_data %>%
  filter(type %in% c("gene", "ncRNA_gene", "miRNA")) %>%
  select(seqnames, start, end, gene_id, strand)

# Adjust start and end based on strand
filtered_gff3_data <- filtered_gff3_data %>%
  mutate(
    start = ifelse(strand == "+", start - 2000, start - 1000),
    end = ifelse(strand == "+", end + 1000, end + 2000)
  )

# Initialize an empty dataframe to store the results
results <- data.frame()

#convert columns to same type of data
filtered_gff3_data$seqnames <- as.character(filtered_gff3_data$seqnames)
combined_df$seqnames <- as.character(combined_df$seqnames)

# Function to find overlaps
find_overlaps <- function(row) {
  gene_ids <- character(0)
  
  matching_genes <- filtered_gff3_data %>%
    filter(seqnames == row$seqnames) %>%
    filter(
      (row$start <= end & row$end >= start) |
        (row$start == (end + 1) & row$end >= start) |
        (row$start <= (end - 1) & row$end == start)
    )
  
  if (nrow(matching_genes) > 0) {
    gene_ids <- paste(matching_genes$gene_id, collapse = ";")
  } else {
    gene_ids <- "none"
  }
  
  result_row <- data.frame(
    seqnames = row$seqnames,
    start = row$start,
    end = row$end,
    gene_id = gene_ids,
    relationship = row$relationship
  )
  return(result_row)
}

# Apply the function to each row in the combined_df data
results <- lapply(1:nrow(combined_df), function(i) find_overlaps(combined_df[i, ]))
results <- do.call(rbind, results)

#Move the gene_ids to the full data frame of peaks
combined_df$bound_genes <- results$gene_id

# Write the results to a CSV file
write.csv(combined_df, "Diffbind Peak with gene Annotation.csv", row.names = FALSE)


###################################################
#Combine the ChIP-seq data with DESeq2 results
#Do not edit this first file. it is included in the folder
genes <- read.csv("gene_ids.csv", header = TRUE, sep = ",")

# Specify the exact relationship for filtering
exact_relationship <- "mutual of samples1&2"

# Filter combined_df based on the exact relationship value
filtered_df <- combined_df %>%
  filter(relationship == exact_relationship)

# Specify the column name for the output in genes dataframe
output_column <- paste(gsub(" ", "_", exact_relationship), "count", sep = "_")

# Function to count occurrences of gene_id in bound_genes and return counts
count_gene_occurrences <- function(gene_id) {
  gene_counts <- str_count(filtered_df$bound_genes, gene_id)
  return(sum(gene_counts))
}

# Apply the function to gene_id and store the counts in the associated column in genes
genes[[output_column]] <- sapply(genes$gene_id, count_gene_occurrences)


#########################
# Specify the exact relationship for filtering
exact_relationship <- "mutual of samples1&2; sample1 high"

# Filter combined_df based on the exact relationship value
filtered_df <- combined_df %>%
  filter(relationship == exact_relationship)

# Specify the column name for the output in genes dataframe
output_column <- paste(gsub(" ", "_", exact_relationship), "count", sep = "_")

# Function to count occurrences of gene_id in bound_genes and return counts
count_gene_occurrences <- function(gene_id) {
  gene_counts <- str_count(filtered_df$bound_genes, gene_id)
  return(sum(gene_counts))
}

# Apply the function to gene_id and store the counts in the associated column in genes
genes[[output_column]] <- sapply(genes$gene_id, count_gene_occurrences)


#########################
# Specify the exact relationship for filtering
exact_relationship <- "mutual of samples1&2; sample2 high"

# Filter combined_df based on the exact relationship value
filtered_df <- combined_df %>%
  filter(relationship == exact_relationship)

# Specify the column name for the output in genes dataframe
output_column <- paste(gsub(" ", "_", exact_relationship), "count", sep = "_")

# Function to count occurrences of gene_id in bound_genes and return counts
count_gene_occurrences <- function(gene_id) {
  gene_counts <- str_count(filtered_df$bound_genes, gene_id)
  return(sum(gene_counts))
}

# Apply the function to gene_id and store the counts in the associated column in genes
genes[[output_column]] <- sapply(genes$gene_id, count_gene_occurrences)


#################################
# Specify the exact relationship for filtering
exact_relationship <- "sample1 unique"

# Filter combined_df based on the exact relationship value
filtered_df <- combined_df %>%
  filter(relationship == exact_relationship)

# Specify the column name for the output in genes dataframe
output_column <- paste(gsub(" ", "_", exact_relationship), "count", sep = "_")

# Function to count occurrences of gene_id in bound_genes and return counts
count_gene_occurrences <- function(gene_id) {
  gene_counts <- str_count(filtered_df$bound_genes, gene_id)
  return(sum(gene_counts))
}

# Apply the function to gene_id and store the counts in the associated column in genes
genes[[output_column]] <- sapply(genes$gene_id, count_gene_occurrences)

########################################
# Specify the exact relationship for filtering
exact_relationship <- "sample2 unique"

# Filter combined_df based on the exact relationship value
filtered_df <- combined_df %>%
  filter(relationship == exact_relationship)

# Specify the column name for the output in genes dataframe
output_column <- paste(gsub(" ", "_", exact_relationship), "count", sep = "_")

# Function to count occurrences of gene_id in bound_genes and return counts
count_gene_occurrences <- function(gene_id) {
  gene_counts <- str_count(filtered_df$bound_genes, gene_id)
  return(sum(gene_counts))
}

# Apply the function to gene_id and store the counts in the associated column in genes
genes[[output_column]] <- sapply(genes$gene_id, count_gene_occurrences)

# Display the updated genes data frame
print(genes)
write.csv(genes, "gene binding events.csv", row.names = FALSE)

###################################
#Begin combination with RNA seq data
#Edit the names of where your RNA seq files are

sample1_rna <- read.csv("./rna/phb_induced_vs_mock.csv", header = TRUE)
sample1_rna <- data.frame("dominant_isoform" = sample1_rna$X, "sample1_log2_fold" = sample1_rna$log2FoldChange, "sample1_padj" = sample1_rna$padj)
sample2_rna <- read.csv("./rna/phb_delta_induced_vs_mock.csv", header = TRUE)
sample2_rna <- data.frame("dominant_isoform" = sample2_rna$X, "sample2_log2_fold" = sample2_rna$log2FoldChange, "sample2_padj" = sample2_rna$padj)

#Merge RNA seq data with gene binding information
ChIP_RNA <- left_join(genes, sample1_rna, by = "dominant_isoform")
ChIP_RNA <- left_join(ChIP_RNA, sample2_rna, by = "dominant_isoform")

#Identify Mututal, sample1, and sample2 target genes
ChIP_RNA$target_of_sample1 <- ifelse(
  rowSums(ChIP_RNA[, 4:7]) > 0 &
    (ChIP_RNA[, 9] >= 1 | ChIP_RNA[, 9] <= -1) &
    ChIP_RNA[, 10] <= 0.1,
  "YES",
  "NO"
)

ChIP_RNA$target_of_sample2 <- ifelse(
  rowSums(ChIP_RNA[, 4:7]) > 0 &
    (ChIP_RNA[, 9] >= 1 | ChIP_RNA[, 9] <= -1) &
    ChIP_RNA[, 10] <= 0.1,
  "YES",
  "NO"
)

ChIP_RNA$target_of_sample2 <- ifelse(
  rowSums(ChIP_RNA[, c(4, 5, 6, 8)]) > 0 &
    (ChIP_RNA[, 11] >= 1 | ChIP_RNA[, 11] <= -1) &
    ChIP_RNA[, 12] <= 0.1,
  "YES",
  "NO"
)

ChIP_RNA$target_of <- ifelse(
  ChIP_RNA[, 13] == "YES" & ChIP_RNA[, 14] == "YES",
  "Both",
  ifelse(
    ChIP_RNA[, 13] == "YES" & ChIP_RNA[, 14] != "YES",
    "Sample1",
    ifelse(
      ChIP_RNA[, 14] == "YES" & ChIP_RNA[, 13] != "YES",
      "Sample2",
      "Neither"
    )
  )
)

#Remove the redundant data
ChIP_RNA <- data.frame(ChIP_RNA[1:12],ChIP_RNA[15])

#Save the ChIP x RNA comparison
write.csv(ChIP_RNA, "true targets.csv", row.names = FALSE)

##############################################
#Generate mutually regulated genes plot
#Mutual Genes
mutual_targets <- subset(ChIP_RNA, target_of == "Both")
g <- ggplot(mutual_targets, aes(x=sample1_log2_fold, y=sample2_log2_fold)) + theme_classic()
g <- g + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_line(linewidth = 0.5, color = "#000000"),
               axis.text.y=element_blank(),
               axis.ticks.y=element_line(linewidth = 0.5, color = "#000000"),
               axis.title =element_blank(),
               axis.line =element_line(linewidth = 0.5),
               axis.ticks.length = unit(0.05, "in"),
               legend.position = "none",
               panel.background = element_rect(fill='transparent'),
               plot.background = element_rect(fill='transparent', color=NA),
)
#Edit the color below to change the color of the points in this plot
g <- g + geom_point(size = 0.1, color = "#FF3399")
#Set limits of the plot axis, this limits it to 8 and squishes any points beyond this point to still be included
g <- g + 
  scale_x_continuous(
    limits = c(-8, 8), 
    oob = scales::squish,
    breaks = seq(-8, 8, by = 2)  # Adjust the breaks to have a major tick every 2 units
  )
g <- g + 
  scale_y_continuous(
    limits = c(-8, 8), 
    oob = scales::squish,
    breaks = seq(-8, 8, by = 2)  # Adjust the breaks to have a major tick every 2 units
  )
g <- g + geom_hline(yintercept = 0)
g <- g + geom_vline(xintercept = 0)
g <- g + geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", size = 0.5)
#If you want to label a specific gene, include its gene_id as the AT number here. Example below labels ZPR4
#g <- g + geom_text_repel(data=mutual_targets[which(mutual_targets$symbol=="AT2G36307"),], aes(label= "ZPR4"), force = 1, min.segment.length = 0, size = 4, segment.size = 1, nudge_x = -3, nudge_y = -3)
g <- g + theme(text=(element_text(size=10, colour = "black", face = "bold")))
g <- g + theme(axis.text = (element_text(colour = "black")))
g
# Save the plot
# X-axis is the log2 fold change of sample1, Y-axis is the log2 fold change of sample 2
ggsave("CNA and PHB mutual genes.pdf", plot = g, height = 2, width = 2, units = "in")

#Get the linear regression information
# Fit linear regression
lm_model <- lm(sample2_log2_fold ~ sample1_log2_fold, data = mutual_targets)

# Extract slope and R-squared
slope <- coef(lm_model)[2]  # Slope
r_squared <- summary(lm_model)$r.squared  # R-squared

# Save slope and R-squared to a text file
output <- data.frame(Number_of_genes = nrow(mutual_targets),Slope = slope, R_squared = r_squared)
write.table(output, "linear_regression_info_of_mutual_targets_plot.txt", sep = "\t", row.names = FALSE)

#############################################
#Generate Correlation Plot
#This plot will filter for genes that are mutually bound, with no specific binding sites
Correlation_data <- ChIP_RNA %>%
  filter((target_of == "Sample1" | target_of == "Sample2") & sample1_unique_count == 0 & sample2_unique_count == 0)
g <- ggplot(Correlation_data, aes(x=sample1_log2_fold, y=sample2_log2_fold, color = target_of)) + theme_classic()
g <- g + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_line(linewidth = 0.5, color = "#000000"),
               axis.text.y=element_blank(),
               axis.ticks.y=element_line(linewidth = 0.5, color = "#000000"),
               axis.title =element_blank(),
               axis.line =element_line(linewidth = 0.5),
               axis.ticks.length = unit(0.05, "in"),
               legend.position = "none",
               panel.background = element_rect(fill='transparent'),
               plot.background = element_rect(fill='transparent', color=NA),
)
g <- g + geom_point(size = 0.1)
g <- g + 
  scale_x_continuous(
    limits = c(-8, 8), 
    oob = scales::squish,
    breaks = seq(-8, 8, by = 2)  # Adjust the breaks to have a major tick every 2 units
  )
g <- g + 
  scale_y_continuous(
    limits = c(-8, 8), 
    oob = scales::squish,
    breaks = seq(-8, 8, by = 2)  # Adjust the breaks to have a major tick every 2 units
  )
g <- g + geom_hline(yintercept = 0)
g <- g + geom_vline(xintercept = 0)
g <- g + theme(text=(element_text(size=10, colour = "black", face = "bold")))
g <- g + theme(axis.text = (element_text(colour = "black")))
# Edit this code below to change the color for the samples
g <- g + scale_color_manual(values = c("Sample1" = "#FF0000", "Sample2" = "#FF66CC"))
g
# Save the plot
# X-axis is the log2 fold change of sample1, Y-axis is the log2 fold change of sample 2
ggsave("CNA and PHB correlation plot.pdf", plot = g, height = 2, width = 2, units = "in")