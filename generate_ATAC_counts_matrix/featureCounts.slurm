#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60000
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
/sfs/qumulo/qproject/CPHG/MILLER/CAD_QTL/coronary_QTL/software/Feature_Count/bin/featureCounts -p -a coronary_consensus_peaks.saf -F SAF -o coronary_artery_counts_matrix \
UVA001_CA_ATAC_sort_dedup.bam \
UVA002_CA_ATAC_sort_dedup.bam \
UVA005_CA_ATAC_sort_dedup.bam \
UVA007_CA_ATAC_sort_dedup.bam \
UVA008_CA_ATAC_2_sort_dedup.bam \
UVA011_CA_ATAC_sort_dedup.bam \
UVA012_CA_ATAC_sort_dedup.bam \
UVA014_CA_ATAC_2_sort_dedup.bam \
UVA015_CA_ATAC_2_sort_dedup.bam \
UVA018_CA_ATAC_1_sort_dedup.bam \
UVA019_CA_ATAC_3_sort_dedup.bam \
UVA020_CA_ATAC_2_sort_dedup.bam \
UVA021_CA_ATAC_sort_dedup.bam \
UVA022_CA_ATAC_2_sort_dedup.bam \
UVA023_CA_ATAC_sort_dedup.bam \
UVA025_CA_ATAC_sort_dedup.bam \
UVA026_CA_ATAC_sort_dedup.bam \
UVA027_CA_ATAC_sort_dedup.bam \
UVA028_CA_ATAC_sort_dedup.bam \
UVA030_CA_ATAC_sort_dedup.bam \
UVA031_CA_ATAC_sort_dedup.bam \
UVA032_CA_ATAC_sort_dedup.bam \
UVA033_CA_ATAC_sort_dedup.bam \
UVA036_CA_ATAC_sort_dedup.bam \
UVA037_CA_ATAC_sort_dedup.bam \
UVA038_CA_ATAC_2_sort_dedup.bam \
UVA039_CA_ATAC_sort_dedup.bam \
UVA040_CA_ATAC_sort_dedup.bam \
UVA041_CA_ATAC_sort_dedup.bam \
UVA044_CA_ATAC_sort_dedup.bam \
UVA045_CA_ATAC_sort_dedup.bam \
UVA046_CA_ATAC_sort_dedup.bam \
UVA047_CA_ATAC_sort_dedup.bam \
UVA049_CA_ATAC_sort_dedup.bam \
UVA050_CA_ATAC_sort_dedup.bam \
UVA051_CA_ATAC_sort_dedup.bam \
UVA052_CA_ATAC_2_sort_dedup.bam \
UVA053_CA_ATAC_sort_dedup.bam \
UVA054_CA_ATAC_sort_dedup.bam \
UVA056_CA_ATAC_sort_dedup.bam \
UVA057_CA_ATAC_3_sort_dedup.bam \
UVA058_CA_ATAC_1_sort_dedup.bam \
UVA060_CA_ATAC_sort_dedup.bam \
UVA061_CA_ATAC_sort_dedup.bam \
UVA062_CA_ATAC_sort_dedup.bam \
UVA063_CA_ATAC_sort_dedup.bam \
UVA064_CA_ATAC_sort_dedup.bam \
UVA065_CA_ATAC_sort_dedup.bam \
UVA066_CA_ATAC_sort_dedup.bam \
UVA067_CA_ATAC_sort_dedup.bam \
UVA068_CA_ATAC_sort_dedup.bam \
UVA069_CA_ATAC_sort_dedup.bam \
UVA070_CA_ATAC_sort_dedup.bam \
UVA071_CA_ATAC_sort_dedup.bam \
UVA072_CA_ATAC_sort_dedup.bam \
UVA073_CA_ATAC_sort_dedup.bam \
UVA074_CA_ATAC_sort_dedup.bam \
UVA078_CA_ATAC_sort_dedup.bam \
UVA079_CA_ATAC_sort_dedup.bam \
UVA080_CA_ATAC_sort_dedup.bam \
UVA081_CA_ATAC_sort_dedup.bam \
UVA082_CA_ATAC_sort_dedup.bam \
UVA084_CA_ATAC_sort_dedup.bam \
UVA085_CA_ATAC_sort_dedup.bam \
UVA086_CA_ATAC_sort_dedup.bam \
UVA087_CA_ATAC_sort_dedup.bam \
UVA088_CA_ATAC_sort_dedup.bam \
UVA091_CA_ATAC_sort_dedup.bam \
UVA094_CA_ATAC_sort_dedup.bam \
UVA096_CA_ATAC_sort_dedup.bam \
UVA098_CA_ATAC_sort_dedup.bam \
UVA099_CA_ATAC_sort_dedup.bam \
UVA100_CA_ATAC_sort_dedup.bam \
UVA101_CA_ATAC_sort_dedup.bam \
UVA102_CA_ATAC_sort_dedup.bam \
UVA103_CA_ATAC_sort_dedup.bam \
UVA104_CA_ATAC_sort_dedup.bam \
UVA106_CA_ATAC_sort_dedup.bam \
UVA107_CA_ATAC_sort_dedup.bam \
UVA108_CA_ATAC_sort_dedup.bam \
UVA111_CA_ATAC_sort_dedup.bam \
UVA112_CA_ATAC_sort_dedup.bam \
UVA113_CA_ATAC_sort_dedup.bam \
UVA115_CA_ATAC_sort_dedup.bam \
UVA116_CA_ATAC_sort_dedup.bam \
UVA119_CA_ATAC_sort_dedup.bam \
UVA120_CA_ATAC_sort_dedup.bam \
UVA121_CA_ATAC_sort_dedup.bam \
UVA122_CA_ATAC_sort_dedup.bam \
UVA123_CA_ATAC_sort_dedup.bam \
UVA129_CA_ATAC_sort_dedup.bam \
UVA130_CA_ATAC_sort_dedup.bam \
UVA131_CA_ATAC_sort_dedup.bam \
UVA132_CA_ATAC_sort_dedup.bam \
UVA133_CA_ATAC_sort_dedup.bam \
UVA134_CA_ATAC_sort_dedup.bam \
UVA136_CA_ATAC_sort_dedup.bam \
UVA138_CA_ATAC_sort_dedup.bam \
UVA140_CA_ATAC_sort_dedup.bam \
UVA141_CA_ATAC_sort_dedup.bam \
UVA146_CA_ATAC_sort_dedup.bam \
UVA150_CA_ATAC_sort_dedup.bam \
UVA151_CA_ATAC_sort_dedup.bam \
UVA152_CA_ATAC_sort_dedup.bam \
UVA155_CA_ATAC_sort_dedup.bam \
UVA156_CA_ATAC_sort_dedup.bam \
UVA157_CA_ATAC_sort_dedup.bam \
UVA159_CA_ATAC_sort_dedup.bam \
UVA161_CA_ATAC_sort_dedup.bam \
UVA163_CA_ATAC_sort_dedup.bam \
UVA164_CA_ATAC_sort_dedup.bam \
UVA165_CA_ATAC_sort_dedup.bam \
UVA167_CA_ATAC_sort_dedup.bam \
UVA168_CA_ATAC_sort_dedup.bam \
UVA170_CA_ATAC_sort_dedup.bam \
UVA171_CA_ATAC_sort_dedup.bam \
UVA173_CA_ATAC_sort_dedup.bam \
UVA174_CA_ATAC_sort_dedup.bam \
UVA175_CA_ATAC_sort_dedup.bam
