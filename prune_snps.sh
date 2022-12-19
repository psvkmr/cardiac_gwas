#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=min_prune
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/min/gwas_prune_%A_%a.out
#SBATCH --array=[1-22]%6


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
#Newer versions throw error about IDs
module load apps/plink2/2.0.0a2

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/min_run
resources=/scratch/users/k2142172/resources

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/sample_filt_chr${i}*


plink2 --pfile ${out_dir}/sample_filt_chr${i} --memory 16000 \
 --indep-pairwise 500 50 0.2 --out ${out_dir}/LDprune_chr${i} 
wait
plink2 --pfile ${out_dir}/sample_filt_chr${i} --memory 16000  \
--exclude bed0 ${resources}/highLDRegionsIDs_GRCh37.txt \
--extract ${out_dir}/LDprune_chr${i}.prune.in --make-pgen --out ${out_dir}/pruned_chr${i}
wait
plink2 --pfile ${out_dir}/pruned_chr${i} --memory 16000 --make-bed --out ${out_dir}/pruned_chr${i}
