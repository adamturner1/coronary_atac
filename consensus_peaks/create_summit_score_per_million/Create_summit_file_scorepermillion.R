#Load packages
library(tidyverse)


#Specify narrowPeak file column names
colnames10 <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")


#Read filtered peak file
#Create summit bed file
#Calculate score per million
UVA001 <- read.delim("UVA001_CA_ATAC_filtered.narrowPeak", header=F, sep="\t", col.names=colnames10) %>% mutate(summit = chromStart + peak) %>% mutate (summit2 = summit) %>% select(chrom, summit, summit2, name, score, strand, signalValue, pValue, qValue, peak) %>% mutate(new_score=pValue/(sum(pValue)/1000000)) %>% select(chrom, summit, summit2, name, new_score, strand)


#Write table
write.table(UVA001, file="UVA001_CA_ATAC_summits_spm.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
