#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=genotype_filt
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/min/genotype_filt.out


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2

resources=/scratch/users/k2142172/resources
geno_data=${resources}/ukb_genotyped/ukb_binary_v2
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/min_run

keep_samples=${out_dir}/sample_filt_chr1.psam

mkdir -p $out_dir

ls -alht ${geno_data}*

#plink2 --bfile ${geno_data} --memory 16000 --threads 2  \
#--keep $keep_samples --make-pgen --out ${out_dir}/genotype_filt

#wait

#plink2 --pfile ${out_dir}/genotype_filt --memory 16000 \
# --indep-pairwise 500 50 0.2 --out ${out_dir}/LDprune_genotype 

#wait 

plink2 --pfile ${out_dir}/genotype_filt --memory 16000 --exclude bed0 ${resources}/highLDRegionsIDs_GRCh37.txt --extract ${out_dir}/LDprune_genotype.prune.in --make-pgen --out ${out_dir}/genotype_pruned 

wait 

plink2 --pfile ${out_dir}/genotype_pruned --memory 16000 --make-bed --out ${out_dir}/genotype_pruned
