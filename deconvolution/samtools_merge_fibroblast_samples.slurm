#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80000
#SBATCH --time=3-00:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load samtools
samtools merge Fibroblast_all_reads.bam Fibroblast_all_bam/UVA005_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA007_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA014_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA016_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA017_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA020_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA024_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA030_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA036_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA039_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA040_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA043_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA046_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA049_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA053_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA054_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA057_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA058_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA064_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA066_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA079_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA080_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA083_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA094_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA108_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA112_Fibroblast_filter.sorted.bam \
Fibroblast_all_bam/UVA124_Fibroblast_filter.sorted.bam
