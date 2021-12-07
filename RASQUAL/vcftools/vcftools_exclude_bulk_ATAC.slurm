#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40000
#SBATCH --time=12:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load vcftools
vcftools --gzvcf UVA_coronary_hg38_merged_imputed.vcf.gz --remove excludedsamples_bulk_ATAC.txt --recode --out UVA_hg38_coronary_bulk_ATAC
