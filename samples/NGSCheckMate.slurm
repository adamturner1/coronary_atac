#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=60000
#SBATCH --time=20:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load anaconda/2019.10-py2.7
module load samtools
module load bcftools
python2 ncm.py -B -d /scratch/amt2ug/bam_sort_dedup -bed SNP/SNP_GRCh38_hg38_wChr.bed -O /scratch/amt2ug/NCM/NGSCheckMate_ATAC
