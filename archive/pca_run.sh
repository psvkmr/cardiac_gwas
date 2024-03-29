#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=gwas_pca
#SBATCH --time=04:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/gwas_pca_%A_%a.out

# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run

mkdir -p $out_dir

ls -alht ${out_dir}/all_chr_snpfilt*

plink2 --pfile ${out_dir}/all_chr_snpfilt --memory 40000 --threads 2 --pca --out ${out_dir}/all_chr_pca
