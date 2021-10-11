#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=80000
#SBATCH --time=4-0:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
module load gcc/9.2.0
module load gsl/2.4
module load rasqual
module load tabix
module load gparallel/20170822
n_lines=$(wc -l bulk_input/rasqual.bulk.input.txt | cut -d" " -f1)
parallel -j10 bash rasqual_bulk.sh \
bulk_input/rasqual.bulk.input.txt {} \
bulk_input/Y.bin \
bulk_input/K.bin \
bulk_input/X.bin \
asvcf/bulk_coronary_ATAC_AS.vcf.gz \
bulk_output ::: `seq $n_lines`
