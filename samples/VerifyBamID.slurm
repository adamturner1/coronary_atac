#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40000
#SBATCH --time=1-12:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load bioconda
verifyBamID --vcf UVA_coronary_hg38_merged_imputed_chr22.vcf.gz --bam UVA001_CA_ATAC_sort_dedup.bam --out UVA001 --verbose --ignoreRG --best
