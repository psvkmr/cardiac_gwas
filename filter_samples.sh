#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=min_samplefilt
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/min/gwas_samplefilt_%A_%a.out
#SBATCH --array=[1-22]%6


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2

remove_samples=$1

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/min_run

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/snp_filt_chr${i}*

plink2 --pfile ${out_dir}/snp_filt_chr${i} --memory 24000 --threads 2  \
--remove $remove_samples --make-pgen --out ${out_dir}/sample_filt_chr${i}

