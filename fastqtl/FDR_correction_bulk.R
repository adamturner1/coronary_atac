#Load packages
library(tidyverse)

#Read results table
d = read.table("permutations.adaptive_allchr.txt", hea=F, stringsAsFactors=F)


#Assign column names
#Newer version of FastQTL has 11 (not 10) columns
colnames(d) = c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")


#Benjamini & Hochberg correction
d$bh = p.adjust(d$bpval, method="fdr")


#Extract all signifcant MP-QTL pairs
write.table(d[which(d$bh <= 0.05), c(1,6)], "permutations.atac.benjamini.fdr05.txt", quote=F, row.names=F, col.names=T, sep="\t")

#write.table(d[which(d$bh <= 0.10), c(1,6)], "permutations.atac.benjamini.fdr10.txt", quote=F, row.names=F, col.names=T, sep="\t")