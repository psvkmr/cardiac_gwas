#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=dpca
#SBATCH --time=4:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/dis/pca_all.out


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/dis_run

mkdir -p $out_dir

ls -alht ${out_dir}/all_pruned_merged*

plink2 --pfile ${out_dir}/all_pruned_merged --memory 40000 --threads 4 \
--pca 20 --out ${out_dir}/pca
