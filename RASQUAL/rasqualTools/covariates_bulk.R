#Load packages
library(tidyverse)

#Bulk covariates
#Covariates are age (binned), sex, and first 3 principal components of ancestry

bulk_covar <- read_csv("bulk_rasqual_covariates_2.csv", col_names = FALSE)
bulk_covar$X1 <- as.factor(bulk_covar$X1)
bulk_covar$X2 <- as.factor(bulk_covar$X2)

#Export covariate table ("X.txt")
#write.table(bulk_covar, file = "X.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
