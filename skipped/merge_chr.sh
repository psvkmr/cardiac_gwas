#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=gwas_merge
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/gwas_merge_%A_%a.out
#SBATCH --array=1-2


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
#module load apps/plink2/2.0.0a2

plink2=/scratch/users/k2142172/packages/plink2
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/pruned_chr${i}*

# breaks if there is another . in the name
ls ${out_dir}/pruned_chr*.pgen | cut -d. -f 1 > ${out_dir}/chrs_to_merge.txt

$plink2 --pmerge-list ${out_dir}/chrs_to_merge.txt pfile --memory 24000 --threads 2  \
--make-pgen --out ${out_dir}/all_chr_merged
