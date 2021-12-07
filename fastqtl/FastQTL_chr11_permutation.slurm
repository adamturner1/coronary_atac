#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60000
#SBATCH --time=2-00:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load gcc/9.2.0
module load fastqtl
fastQTL --vcf genotypes.vcf.gz  --bed  phenotypes.bed.gz  --region  chr11   --window 5e4  --permute 1000 100000  --out permutations.adaptive_chr11.txt.gz  --cov covariates.txt.gz
