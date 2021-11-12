#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=40000
#SBATCH --time=2:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
# HOMER is already loaded
annotatePeaks.pl coronary_artery_consensus_peaks.bed hg38 > coronary_artery_consensus_peaks_annotated.txt
