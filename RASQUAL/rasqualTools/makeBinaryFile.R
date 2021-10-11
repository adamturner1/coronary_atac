#Load packages
library(tidyverse)

#Follow the RASQUAL script
xtxt = "X.txt"
if(!is.na(xtxt)){x=read.table(xtxt,as.is=T)}
if(!is.na(xtxt)){xbin=gsub("txt", "bin", xtxt)}
if(!is.na(xtxt)){fxbin=file(xbin,"wb")}
if(!is.na(xtxt)){writeBin(as.double(c(as.matrix(x))), fxbin)}
if(!is.na(xtxt)){close(fxbin)}
