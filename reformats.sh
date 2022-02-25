#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=gwas_reformat
#SBATCH --time=04:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/gwas_reformat_%A_%a.out
#SBATCH --array=[1-22]%4


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run
bgzip=/scratch/users/k2142172/packages/anaconda3/envs/peddy/bin/bgzip
bcftools=/scratch/users/k2142172/packages/bcftools-1.14/bcftools

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/sample_filt_chr${i}*

plink2 --pfile ${out_dir}/sample_filt_chr${i} --make-bed --out ${out_dir}/sample_filt_chr${i}
plink2 --pfile ${out_dir}/sample_filt_chr${i} --export vcf vcf-dosage=DS-force --out ${out_dir}/sample_filt_chr${i}
wait
# compress the VCF file version of final data files
$bgzip -i ${out_dir}/sample_filt_chr${i}.vcf
#wait
# index created compressed VCF file
#$bcftools index ${out_dir}/sample_filt_chr${i}.vcf.gz

