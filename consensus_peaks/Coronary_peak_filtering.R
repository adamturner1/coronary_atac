#Load packages
library(tidyverse)
library(data.table)


#Read ATAC peak matrices (from intersections)
peak_files <- list.files(pattern="*intersected.narrowPeak")

all_peaks <- lapply(peak_files, fread)


#Use lapply to filter based on row sums (>=2). This step keeps the peaks that are reproducible in other samples
reproducible_peaks <- lapply(all_peaks, function(x) x[rowSums(x[,11:128]) >=2, 1:10])


#Name each item in the list
names(reproducible_peaks) <- c("UVA001", "UVA002", "UVA005", "UVA007", "UVA008", "UVA011", "UVA012", "UVA014", "UVA015", "UVA018", "UVA019", "UVA020", "UVA021", "UVA022", "UVA023", "UVA025", "UVA026", "UVA027", "UVA028", "UVA030", "UVA031", "UVA032", "UVA033", "UVA036", "UVA037", "UVA038", "UVA039", "UVA040", "UVA041", "UVA044", "UVA045", "UVA046", "UVA047", "UVA049", "UVA050", "UVA051", "UVA052", "UVA053", "UVA054", "UVA056", "UVA057", "UVA058", "UVA060", "UVA061", "UVA062", "UVA063", "UVA064", "UVA065", "UVA066", "UVA067", "UVA068", "UVA069", "UVA070", "UVA071", "UVA072", "UVA073", "UVA074", "UVA078", "UVA079", "UVA080", "UVA081", "UVA082", "UVA084", "UVA085", "UVA086", "UVA087", "UVA088", "UVA091", "UVA094", "UVA096", "UVA098", "UVA099", "UVA100", "UVA101", "UVA102", "UVA103", "UVA104", "UVA106", "UVA107", "UVA108", "UVA111", "UVA112", "UVA113", "UVA115", "UVA116", "UVA119", "UVA120", "UVA121", "UVA122", "UVA123", "UVA129", "UVA130", "UVA131", "UVA132", "UVA133", "UVA134", "UVA136", "UVA138", "UVA140", "UVA141", "UVA146", "UVA150", "UVA151", "UVA152", "UVA155", "UVA156", "UVA157", "UVA159", "UVA161", "UVA163", "UVA164", "UVA165", "UVA167", "UVA168", "UVA170", "UVA171", "UVA173", "UVA174", "UVA175")


#Write filtered peaks to new files
lapply(names(reproducible_peaks), function(x) write.table(reproducible_peaks[[x]], file=paste0(x,'_CA_ATAC_filtered.narrowPeak'), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE))
