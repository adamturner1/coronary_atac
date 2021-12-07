#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=60000
#SBATCH --time=7-0:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab-VIP
module load gcc/9.2.0
module load samtools
module load tabix
module load rasqual
bash ./createASVCF.sh paired_end bam.list.bulk.txt UVA_hg38_coronary_bulk_ATAC_filtered.vcf.gz bulk_coronary_ATAC_AS.vcf.gz atac
