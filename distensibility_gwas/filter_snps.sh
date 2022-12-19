#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=dis_snpfilt
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/dis/gwas_snpfilt_%A_%a.out
#SBATCH --array=[1-22]%6


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/dis_run
resources=/scratch/users/k2142172/resources

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/bgen_filt_chr${i}*

plink2 --pfile ${out_dir}/bgen_filt_chr${i} --memory 24000 --threads 2  \
--maf 0.01 --geno 0.05 --hwe 1e-6 --snps-only --make-pgen --extract ${resources}/rsIDs_In_HRC_Or_Genotyped.txt --out ${out_dir}/snp_filt_chr${i}

wait

plink2 --pfile ${out_dir}/snp_filt_chr${i} --memory 24000 --threads 2 \
 --freq counts --out ${out_dir}/snp_filt_chr${i}_freq 
