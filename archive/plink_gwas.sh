#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=gwas_plink
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/gwas_plink_%A_%a.out
#SBATCH --array=1-2


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

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/qc3_chr${i}_snpfilt*

plink2 --pfile ${out_dir}/qc3_chr${i}_snpfilt \
 --memory 24000 \
--threads 1  \
--pheno /scratch/users/k2142172/outputs/cardiac_gwas/res_distensibility3.txt \
--covar ${out_dir}/qc_chr${i}_pca.eigenvec \
--pheno-name res_distensibility \
--covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--glm hide-covar \
--ci 0.95 \
--adjust-file \
--out ${out_dir}/gwas_results_chr${i}
