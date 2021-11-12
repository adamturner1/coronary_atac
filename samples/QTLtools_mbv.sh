#!/bin/bash

#Load modules
module load gcc/9.2.0
module load qtltools

cat samples.txt | while read sample || [[ -n $line ]];
do
        sbatch -A cphg-millerlab-VIP \
               -p parallel \
               -t 12:00:00 \
               -N 2 \
               -n 20 \
               --mem=60000 \
               --wrap "QTLtools mbv --bam /scratch/amt2ug/bam_sort_dedup/${sample} \
                                    --vcf /scratch/amt2ug/ASEReadCounter/UVA_coronary_hg38_merged_imputed.vcf.gz \
                                    --out match_out_${sample}"
done
