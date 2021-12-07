#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60000
#SBATCH --time=3-00:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load gcc/9.2.0
module load fastqtl
fastQTL --vcf genotypes.vcf.gz  --bed  phenotypes.bed.gz  --region  chr5   --window 5e4   --out nominals.default_chr5.txt.gz  --cov covariates.txt.gz
