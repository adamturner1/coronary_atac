#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80000
#SBATCH --time=24:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load samtools
samtools index Endothelial_all_reads.bam
samtools index Fibroblast_all_reads.bam
samtools index Macrophage_all_reads.bam
samtools index Mast_all_reads.bam
samtools index Pericyte_all_reads.bam
samtools index Plasma_all_reads.bam
samtools index SMC_all_reads.bam
samtools index T_all_reads.bam
samtools index unknown_all_reads.bam
