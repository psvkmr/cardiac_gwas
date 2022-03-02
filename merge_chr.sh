#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=gwas_merge
#SBATCH --time=00:30:00
#SBATCH --mem=5G
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/merge.out


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
#module load apps/plink2/2.0.0a2

plink2=/scratch/users/k2142172/packages/plink2
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/clean_gwas

mkdir -p $out_dir

ls ${out_dir}/pruned_chr*.pgen | cut -d. -f 1 > ${out_dir}/pruned_to_merge.txt
wait
$plink2 --pmerge-list ${out_dir}/pruned_to_merge.txt pfile --memory 5000 \
--merge-max-allele-ct 2 --out ${out_dir}/all_pruned_merged
