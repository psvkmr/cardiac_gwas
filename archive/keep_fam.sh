#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=gwas_keepfam
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/gwas_keepfam_%A_%a.out
#SBATCH --array=1-2


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2

raw_data=/scratch/users/k2142172/outputs/cardiac_gwas/ukb_files
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/converted_chr${i}*


 plink2 --pfile ${out_dir}/converted_chr${i} --memory 24000 --threads 1 --keep-fam ${raw_data}/IDs_final.txt --make-pgen \
--out ${out_dir}/qc1_chr${i}_keepfam
