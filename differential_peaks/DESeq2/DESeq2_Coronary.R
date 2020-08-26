#Load packages
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)

#Read the counts matrix (cts)
cts <- read.table("coronary_artery_counts_matrix", header=T, row.names=1)
cts <- cts[,6:124]

#Read the metadata (coldata)
coldata <- read_csv("atac_metadata.csv", col_names = TRUE)

coldata$Classification <- as.factor(coldata$Classification)
coldata$Sex <- as.factor(coldata$Sex)

#Shorten names
sample_names <- coldata$Name
rownames(coldata) <- sample_names
colnames(cts) <- sample_names

#Perform DESeq
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~FRiP + Sex + Classification)
dds <- DESeq(dds)

#Get results
resultsNames(dds)
res <- results(dds, name="Classification_Normal_vs_Ischemic")

#Generate MA plot
plotMA(res, ylim=c(-2,2))

#Visualize using PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Classification"))

#Generate volcano plot
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue', FCcutoff = 0.58)
