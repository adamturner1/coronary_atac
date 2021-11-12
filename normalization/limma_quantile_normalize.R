#Load packages
library(tidyverse)
library(limma)

#Read raw counts
a <- read.table("bulk.counts2.txt", header=F, row.names=1)


#Normalize columns to have the same quantiles (limma package)
bulk_norm <- normalizeQuantiles(a)


#Assign column names
colnames(bulk_norm) <- c("UVA001", "UVA002", "UVA005", "UVA007", "UVA008", "UVA011", "UVA012", "UVA014",
                    "UVA015", "UVA018", "UVA019", "UVA020", "UVA021", "UVA022", "UVA023", "UVA025",
                    "UVA026", "UVA027", "UVA028", "UVA030", "UVA031", "UVA032", "UVA033", "UVA036",
                    "UVA037", "UVA038", "UVA039", "UVA040", "UVA041", "UVA044", "UVA045", "UVA046",
                    "UVA047", "UVA049", "UVA050", "UVA051", "UVA052", "UVA053", "UVA054", "UVA056",
                    "UVA057", "UVA058", "UVA059", "UVA060", "UVA061", "UVA062", "UVA063", "UVA064",
                    "UVA065", "UVA066", "UVA067", "UVA068", "UVA069", "UVA070", "UVA071", "UVA072",
                    "UVA073", "UVA074", "UVA078", "UVA079", "UVA080", "UVA081", "UVA082", "UVA084",
                    "UVA085", "UVA086", "UVA087", "UVA088", "UVA091", "UVA094", "UVA096", "UVA098",
                    "UVA099", "UVA100", "UVA101", "UVA102", "UVA103", "UVA104", "UVA106", "UVA107",
                    "UVA108", "UVA111", "UVA112", "UVA113", "UVA116", "UVA119", "UVA120", "UVA121",
                    "UVA122", "UVA123", "UVA129", "UVA130", "UVA132", "UVA133", "UVA136", "UVA138",
                    "UVA140", "UVA141", "UVA146", "UVA150", "UVA151", "UVA152", "UVA156", "UVA157",
                    "UVA159", "UVA161", "UVA163", "UVA165", "UVA167", "UVA168", "UVA170", "UVA171",
                    "UVA174", "UVA175")

#Round
bulk_norm <- round(bulk_norm, 3)

#Add extra columns
peaks <- rownames(bulk_norm)
bulk_norm <- bulk_norm %>% mutate(peak_id = peaks)
bulk_norm <- bulk_norm %>% mutate(ID = peak_id)
bulk_norm <- bulk_norm %>% separate(peak_id, c("Chr", "start", "end"))

#Change the column order to get in bed format for FastQTL
bulk_norm <- bulk_norm %>% relocate(Chr, .before = UVA001)
bulk_norm <- bulk_norm %>% relocate(start, .after = Chr)
bulk_norm <- bulk_norm %>% relocate(end, .after = start)
bulk_norm <- bulk_norm %>% relocate(ID, .after = end)

#Write table
#write.table(bulk_norm, file="bulk_quantile_normalized_counts_FastQTL.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

#Write table specifically for FastQTL
write.table(bulk_norm, file="bulk_quantile_normalized_counts_FastQTL_phenotypes.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
