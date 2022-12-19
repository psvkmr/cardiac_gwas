#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=dis_qcsample
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/dis/gwas_qcsample_%A_%a.out
#SBATCH --array=[1-22]%6


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2
plink2=/scratch/users/k2142172/packages/plink2
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/dis_run

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/snp_filt_chr${i}*

$plink2 --pfile ${out_dir}/snp_filt_chr${i} --memory 24000 --threads 2 --het --out ${out_dir}/snp_filt_chr${i}_het
$plink2 --pfile ${out_dir}/snp_filt_chr${i} --memory 24000 --threads 2 --missing --out ${out_dir}/snp_filt_chr${i}_miss

