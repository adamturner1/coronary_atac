#Load packages
library(rasqualTools)
library(tidyverse)
library(data.table)

#Import bulk ATAC count matrix
bulk_atac_counts <- fread("coronary_artery_counts_matrix_raw_063021_2")
bulk_atac_counts <- as_tibble(bulk_atac_counts)


#Get GC content of peaks
gc_peaks <- fread("all_peaks_gc.bed")


#Add GC content column
bulk_atac_counts <- cbind(bulk_atac_counts, gc_peaks$`6_pct_gc`)


#Remove pediatric samples
#UVA173 is same as UVA059
colnames(bulk_atac_counts) <- c("Geneid", "Chr", "Start", "End", "Strand", "Length",
                                "UVA001", "UVA002", "UVA005", "UVA007", "UVA008", "UVA011", "UVA012", "UVA014", "UVA015", "UVA018",
                                "UVA019", "UVA020", "UVA021", "UVA022", "UVA023", "UVA025", "UVA026", "UVA027", "UVA028", "UVA030",
                                "UVA031", "UVA032", "UVA033", "UVA036", "UVA037", "UVA038", "UVA039", "UVA040", "UVA041", "UVA044",
                                "UVA045", "UVA046", "UVA047", "UVA049", "UVA050", "UVA051", "UVA052", "UVA053", "UVA054", "UVA056",
                                "UVA057", "UVA058", "UVA060", "UVA061", "UVA062", "UVA063", "UVA064", "UVA065", "UVA066", "UVA067",
                                "UVA068", "UVA069", "UVA070", "UVA071", "UVA072", "UVA073", "UVA074", "UVA078", "UVA079", "UVA080",
                                "UVA081", "UVA082", "UVA084", "UVA085", "UVA086", "UVA087", "UVA088", "UVA091", "UVA094", "UVA096",
                                "UVA098", "UVA099", "UVA100", "UVA101", "UVA102", "UVA103", "UVA104", "UVA106", "UVA107", "UVA108",
                                "UVA111", "UVA112", "UVA113", "UVA115", "UVA116", "UVA119", "UVA120", "UVA121", "UVA122", "UVA123",
                                "UVA129", "UVA130", "UVA131", "UVA132", "UVA133", "UVA134", "UVA136", "UVA138", "UVA140", "UVA141",
                                "UVA146", "UVA150", "UVA151", "UVA152", "UVA155", "UVA156", "UVA157", "UVA159", "UVA161", "UVA163",
                                "UVA164", "UVA165", "UVA167", "UVA168", "UVA170", "UVA171", "UVA059", "UVA174", "UVA175", "GC")

bulk_atac_counts <- bulk_atac_counts %>% relocate(UVA059, .after = UVA058)

bulk_atac_counts <- bulk_atac_counts %>% select(-UVA115, -UVA131, -UVA134, -UVA155, -UVA164)


#Filter counts matrix - mean read count of 10 or more
bulk_atac_counts <- bulk_atac_counts %>% mutate(mean = rowMeans(bulk_atac_counts[,7:120], na.rm=TRUE)) %>% filter(mean >= 10)


#Prepare gene_data (Geneid is the peak name)
gene_data = dplyr::select(bulk_atac_counts, Geneid, Chr, Strand, Start, End, GC)
gene_data_offset <- gene_data %>% select(Geneid, GC)
colnames(gene_data_offset) <- c("gene_id", "percentage_gc_content")
print(gene_data_offset)


#Format counts matrix
peaknames <- bulk_atac_counts$Geneid
bulk_atac_counts2 <- bulk_atac_counts[, 7:120]
row.names(bulk_atac_counts2) <- peaknames


#Save counts matrix
saveRasqualMatrices(list(bulk = bulk_atac_counts2), "/scratch/amt2ug/rasqualTools_bulk", file_suffix = "counts")


#Calculate size factors
size_factors = rasqualCalculateSampleOffsets(bulk_atac_counts2, gene_data_offset, gc_correct = TRUE)
saveRasqualMatrices(list(bulk = size_factors), "/scratch/amt2ug/rasqualTools_bulk", file_suffix = "size_factors")


#Calculate the number of SNPs overlapping each peak
snp_coords <- read.table("snp_list_4.txt", header = TRUE, sep = "\t")
colnames(snp_coords) <- c("chr", "pos", "snp_id")

gene_metadata <- gene_data[ , 1:5]
colnames(gene_metadata) <- c("gene_id", "chr", "strand", "start", "end")

snp_counts = countSnpsOverlapingPeaks(gene_metadata, snp_coords, cis_window = 1e4)


#Format columns for RASQUAL input
snp_counts <- snp_counts %>% unite("z", chromosome_name, range_start, sep = ":", remove = FALSE) %>% unite("region", z, range_end, sep = "-", remove = FALSE)
snp_counts <- snp_counts %>% mutate(gene_name = gene_id)


#Change order of columns
#Exon start and ends correspond to peak start and end positions
snp_counts_2 <- snp_counts %>% select(gene_id, gene_name, region, cis_snp_count, feature_snp_count, exon_starts, exon_ends)


#Export table
write.table(snp_counts_2, file = "rasqual.bulk.input.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
