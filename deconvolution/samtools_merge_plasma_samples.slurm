#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80000
#SBATCH --time=1-00:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load samtools
samtools merge Plasma_all_reads.bam Plasma/scA_Plasma_filter.sorted.bam \
Plasma/scB_Plasma_filter.sorted.bam \
Plasma/scC_Plasma_filter.sorted.bam \
Plasma/scD_Plasma_filter.sorted.bam \
Plasma/scE_Plasma_filter.sorted.bam \
Plasma/scF_Plasma_filter.sorted.bam \
Plasma/scG_Plasma_filter.sorted.bam \
Plasma/scH_Plasma_filter.sorted.bam \
Plasma/scI_Plasma_filter.sorted.bam \
Plasma/scJ_Plasma_filter.sorted.bam \
Plasma/scK_Plasma_filter.sorted.bam \
Plasma/scL_Plasma_filter.sorted.bam \
Plasma/scM_Plasma_filter.sorted.bam \
Plasma/scN_Plasma_filter.sorted.bam \
Plasma/scO_Plasma_filter.sorted.bam \
Plasma/scP_Plasma_filter.sorted.bam \
Plasma/scQ_Plasma_filter.sorted.bam \
Plasma/scR_Plasma_filter.sorted.bam \
Plasma/scS_Plasma_filter.sorted.bam \
Plasma/scT_Plasma_filter.sorted.bam \
Plasma/scU_Plasma_filter.sorted.bam \
Plasma/scV_Plasma_filter.sorted.bam \
Plasma/scW_Plasma_filter.sorted.bam \
Plasma/scX_Plasma_filter.sorted.bam \
Plasma/scY_Plasma_filter.sorted.bam \
Plasma/scZ_Plasma_filter.sorted.bam \
Plasma/scAA_Plasma_filter.sorted.bam \
Plasma/scAB_Plasma_filter.sorted.bam \
Plasma/scAC_Plasma_filter.sorted.bam \
Plasma/scAD_Plasma_filter.sorted.bam \
Plasma/scAE_Plasma_filter.sorted.bam \
Plasma/scAF_Plasma_filter.sorted.bam \
Plasma/scAG_Plasma_filter.sorted.bam \
Plasma/scAH_Plasma_filter.sorted.bam \
Plasma/scAI_Plasma_filter.sorted.bam \
Plasma/scAJ_Plasma_filter.sorted.bam \
Plasma/scAK_Plasma_filter.sorted.bam \
Plasma/scAL_Plasma_filter.sorted.bam \
Plasma/scAM_Plasma_filter.sorted.bam \
Plasma/scAN_Plasma_filter.sorted.bam \
Plasma/scAO_Plasma_filter.sorted.bam \
Plasma/scAP_Plasma_filter.sorted.bam \
Plasma/scAQ_Plasma_filter.sorted.bam \
Plasma/scAR_Plasma_filter.sorted.bam
