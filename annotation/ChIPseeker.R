#Load packages
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#Read peak file
peak <- readPeakFile("coronary_artery_consensus_peaks.bed")

#Annotate peaks
peakAnno <- annotatePeak(peak, tssRegion = c(-1000, 100), TxDb = txdb, annoDb = "org.Hs.eg.db")

#Plot pie chart
plotAnnoPie(peakAnno)
