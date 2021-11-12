#This script is to export promoter and non-promoter peaks using ChIPseeker
library(ChIPseeker)
library(tidyverse)


#Convert peak annotation to data frame
peak2 <- as.data.frame(peakAnno)


#Get promoter and non-promoter regions
promoter_peaks <- peak2[peak2$annotation == "Promoter",]
no_promoter_peaks <- peak2[peak2$annotation != "Promoter",]
