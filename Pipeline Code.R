#Install the packages below if you dont already have them. Some may be through bioconductor
#Library Packages
library(DiffBind)
library(profileplyr)
library(dplyr)
library(ggplot2)
library(reshape)
library(stringr)
library(rtracklayer)

#Read file that specifies file locations and information for DiffBind
#To keep color schemes consistent later, highly recommend putting your mutant / comparing construct 2nd in the data sheet, and wt as the 3rd
#This should be the only file you need to make changes, change sample names in the file and paths to the appropriate files
#then changing the file path below to where the "files.csv" is saved to
data_files <- read.csv("./outputs/cna v delta v wt/files.csv", sep = ",", header = TRUE)

#Perform DiffBind Analysis
#Common error seen with dba.count: "Read block operation failed with error 1 after 0 of 4 bytes"
#Online says this error occurs when the bam file is corrupted during download
#IDK why, just keep redownloading it until it starts working for this file
#Has taken 10 tries before with me, but if it really is not working for a control, switch it to a control input that has worked
#This code is set up for an analysis of 2 genotypes, with 3 replicates each. If your analysis is different from that, you will need to adjust some of the code 
Comparison <- dba(sampleSheet= data_files)
Comparison <- dba.count(Comparison, summits = 200, minOverlap = 2)
Comparison <- dba.normalize(Comparison)
Comparison <- dba.contrast(Comparison, minMembers = 2)
Comparison <- dba.analyze(Comparison)

#Generate the report of the DiffBind analysis, convert to data frame, arrange in same order
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


#Combine information from the data sets
db <- data.frame(ChIP.DB.1[1:5], ChIP.DB.1[7:8], ChIP.DB.2[8],ChIP.DB.1[11], ChIP.DB.2[11],ChIP.DB.3[11], ChIP.DB.1[12:13], ChIP.DB.2[13], ChIP.DB.1[14])
db1 <- db 

#Identify sites we call as true binding, with at least 2 replicates in one genotype having the same called peak
db <- filter(db, db[12]>=2 | db[13]>=2 | db[14]>=2)

#Begin identifying relationships
#Dont be afraid of errors, they very well may happen if there are no matches. just keep going down
#Identify mutual binding sites
mutual <- filter(db, db[12]>=1 & db[13]>=1 & db[14]>=1 & db[9]>=0.1 & db[10]>=0.1 & db[11]>=0.1)
mutual$relationship <- "true_mutual"

#Identify mutual with higher affinity in one
mutual_1_high <- filter(db, db[12]>=1 & db[13]>=1 & db[14]>=1 & db[9] <= 0.05 & db[10]<=0.05 & db[6] > db[7] & db[6] > db[8])
mutual_1_high$relationship <- "mutual sample1 high"

mutual_2_high <- filter(db, db[12]>=1 & db[13]>=1 & db[14]>=1 & db[9] <= 0.05 & db[11]<=0.05 & db[7] > db[6] & db[7] > db[8])
mutual_2_high$relationship <- "mutual sample2 high"

mutual_3_high <- filter(db, db[12]>=1 & db[13]>=1 & db[14]>=1 & db[10] <= 0.05 & db[11]<=0.05 & db[8] > db[6] & db[8] > db[7])
mutual_3_high$relationship <- "mutual sample3 high"

#Identify mutual with higher affinity in two
mutual_1_2_high <- filter(db, db[12]>=1 & db[13]>=1 & db[14]>=1 & db[9] >= 0.1 & db[10]<=0.05 & db[11] <= 0.05 & db[6] > db[8] & db[7] > db[8])
mutual_1_2_high$relationship <- "mutual samples1&2 high"

mutual_1_3_high <- filter(db, db[12]>=1 & db[13]>=1 & db[14]>=1 & db[9] <= 0.05 & db[10]>=0.1 & db[11] <= 0.05 & db[6] > db[7] & db[8] > db[7])
mutual_1_3_high$relationship <- "mutual samples1&3 high"

mutual_2_3_high <- filter(db, db[12]>=1 & db[13]>=1 & db[14]>=1 & db[9] <= 0.05 & db[10]<=0.05 & db[11] >= 0.1 & db[7] > db[6] & db[8] > db[6])
mutual_2_3_high$relationship <- "mutual samples2&3 high"

#Identify unique binding sites
sample1_unique <- filter(db, db[12]>=2 & db[13] == 0 & db[14]==0 & db[9]<= 0.05, db[10]<=0.05 & db[6] > db[7] & db[6] > db[8])
sample1_unique$relationship <- "sample1 unique"

sample2_unique <- filter(db, db[12]==0 & db[13] >= 2 & db[14]==0 & db[9]<= 0.05, db[11]<=0.05 & db[7] > db[6] & db[7] > db[8])
sample2_unique$relationship <- "sample2 unique"

sample3_unique <- filter(db, db[12]==0 & db[13] == 0 & db[14]>=2 & db[10]<= 0.05, db[11]<=0.05 & db[8] > db[7] & db[8] > db[7])
sample3_unique$relationship <- "sample3 unique"

#Identify mutual binding sites in 2 of the 3 samples (3)
mutual_of_1_2 <- filter(db, db[12]>=1 & db[13]>=1 & db[14]==0 & db[9]>= 0.1 & db[10]<= 0.05 & db[11]<= 0.05 & db[6] > db[8] & db[7] > db[8])
mutual_of_1_2$relationship <- "mutual of samples1&2"

mutual_of_1_3 <- filter(db, db[12]>=1 & db[13]==0 & db[14]>=1 & db[9]<= 0.05 & db[10]>= 0.1 & db[11]<= 0.05 & db[6] > db[7] & db[8] > db[7])
mutual_of_1_3$relationship <- "mutual of samples1&3"

mutual_of_2_3 <- filter(db, db[12]==0 & db[13]>=1 & db[14]>=1 & db[9]<= 0.05 & db[10]<= 0.05 & db[11]>= 0.1 & db[7] > db[6] & db[8] > db[6])
mutual_of_2_3$relationship <- "mutual of samples2&3"

#Identify mutual binding sites in 2 of the 3 samples, one higher than other (6)
mutual_of_1_2_1high <- filter(db, db[12]>=1 & db[13]>=1 & db[14]==0 & db[9]<= 0.05 & db[10]<= 0.05 & db[11]<= 0.05 & db[6] > db[7] & db[6] > db[8])
mutual_of_1_2_1high$relationship <- "mutual of samples1&2; sample1 high"

mutual_of_1_2_2high <- filter(db, db[12]>=1 & db[13]>=1 & db[14]==0 & db[9]<= 0.05 & db[10]<= 0.05 & db[11]<= 0.05 & db[7] > db[6] & db[7] > db[8])
mutual_of_1_2_2high$relationship <- "mutual of samples1&2; sample2 high"

mutual_of_1_3_1high <- filter(db, db[12]>=1 & db[13]==0 & db[14]>=1 & db[9]<= 0.05 & db[10]<= 0.05 & db[11]<= 0.05 & db[6] > db[7] & db[6] > db[8])
mutual_of_1_3_1high$relationship <- "mutual of samples1&3; sample1 high"

mutual_of_1_3_3high <- filter(db, db[12]>=1 & db[13]==0 & db[14]>=1 & db[9]<= 0.05 & db[10]<= 0.05 & db[11]<= 0.05 & db[8] > db[7] & db[8] > db[6])
mutual_of_1_3_3high$relationship <- "mutual of samples1&3; sample3 high"

mutual_of_2_3_2high <- filter(db, db[12]==0 & db[13]>=1 & db[14]>=1 & db[9]<= 0.05 & db[10]<= 0.05 & db[11]>= 0.05 & db[7] > db[6] & db[7] > db[8])
mutual_of_2_3_2high$relationship <- "mutual of samples2&3; sample2 high"

mutual_of_2_3_3high <- filter(db, db[12]==0 & db[13]>=1 & db[14]>=1 & db[9]<= 0.05 & db[10]<= 0.05 & db[11]>= 0.05 & db[8] > db[6] & db[8] > db[7])
mutual_of_2_3_3high$relationship <- "mutual of samples2&3; sample3 high"

#Identify removed peaks (FDRs, wrong binding site (example only 1 in each), confusing relationships, the works)

#combine the dataframes
df_list <- list(mutual, mutual_1_high, mutual_2_high,mutual_3_high,mutual_1_2_high, mutual_1_3_high, mutual_2_3_high,sample1_unique, sample2_unique, sample3_unique, mutual_of_1_2, mutual_of_1_3, mutual_of_2_3, mutual_of_1_2_1high, mutual_of_1_2_2high, mutual_of_1_3_1high, mutual_of_1_3_3high, mutual_of_2_3_2high, mutual_of_2_3_3high)
combined_df <- do.call(rbind, df_list)


# List the row names of binding sites with asigned relationships
rows_to_remove <- c(combined_df$row_name)

# Filter out rows from the origional dffbind results to identify
# sites with no established relationship (q-values wrong, complex relationships)
# adds "removed" as relationship
removed <- db1 %>%
  filter(!row_name %in% rows_to_remove)
removed$relationship <- "removed"

#add the removed sites to the rest of the data
list(mutual, mutual_1_high)
combined_df <- do.call(rbind, list(combined_df, removed))

#output full labeled dataset and quick summary
combined_df <- combined_df[order(combined_df$row_name), ]
write.csv(combined_df, "diffbind_ouput.csv", row.names = FALSE)
summary_df <- data.frame("mutual" = nrow(mutual),"mutual sample1 high" = nrow(mutual_1_high),"mutual sample2 high" = nrow(mutual_2_high) ,"mutual sample3 high" = nrow(mutual_3_high), "mutual samples1&2 high" = nrow(mutual_1_2_high), "mutual samples1&3 high" = nrow(mutual_1_3_high), "mutual samples2&3 high" = nrow(mutual_2_3_high),"sample1 unique" =  nrow(sample1_unique), "sample2 unique" = nrow(sample2_unique), "sample3 unique" = nrow(sample3_unique), "mutual of samples1&2" = nrow(mutual_of_1_2), "mutual of samples1&3" = nrow(mutual_of_1_3), "mutual of samples2&3" = nrow(mutual_of_2_3), "mutual of samples1&2; sample1 high" = nrow(mutual_of_1_2_1high), "mutual of samples1&2; sample2 high" = nrow(mutual_of_1_2_2high),"mutual of samples1&3; sample1 high" = nrow(mutual_of_1_3_1high), "mutual of samples1&3; sample3 high" = nrow(mutual_of_1_3_3high), "mutual of samples2&3; sample2 high" =  nrow(mutual_of_2_3_2high), "mutual of samples2&3; sample3 high" = nrow(mutual_of_2_3_3high), "removed" = nrow(removed))
write.csv(summary_df, "relationship_summary.csv", row.names = FALSE)


###########################################################
#Begin to Generate Heatmaps
#This code assumes that "sample3" is wildtype
#There are 2 codes that can alter the color (scale_fill_gradiant & scale_fill_distiller)
#Distiller will make the scale itself for the most part except the NA color
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
#IDK what im doing from here. I want to start combination with the RNAseq
# The next like 100 lines are sloppy and dont loop but work
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

sample1_rna <- read.csv("./rna/cna_induced_vs_mock.csv", header = TRUE)
sample1_rna <- data.frame("dominant_isoform" = sample1_rna$X, "sample1_log2_fold" = sample1_rna$log2FoldChange, "sample1_padj" = sample1_rna$padj)
sample2_rna <- read.csv("./rna/cna_delta_induced_vs_mock.csv", header = TRUE)
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