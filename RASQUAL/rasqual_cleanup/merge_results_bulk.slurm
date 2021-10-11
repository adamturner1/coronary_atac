#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40000
#SBATCH --time=1:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg-millerlab
for i in {1..22}
do
        cat bulk_output/chr${i}.* > Bulk_caqtl_results_chr${i}.txt
done
