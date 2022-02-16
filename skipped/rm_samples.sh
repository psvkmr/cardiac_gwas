#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=gwas_rmsample
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/gwas_rmsample_%A_%a.out
#SBATCH --array=1-2


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2
plink=/scratch/users/k2142172/packages/anaconda3/bin/plink
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run
env=/scratch/users/k2142172/packages/anaconda/envs/peddy

#conda activate $env
conda activate peddy

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/bgen_filt_chr${i}*

#$plink --bfile ${out_dir}/bgen_filt_chr${i} --memory 24000 --threads 2 --check-sex --out ${out_dir}/bgen_filt_chr${i}_sex_check
#wait
#$plink --bfile ${out_dir}/bgen_filt_chr${i} --memory 24000 --threads 2 --het --out ${out_dir}/bgen_filt_chr${i}_het_check
#wait
#$plink --bfile ${out_dir}/bgen_filt_chr${i} --memory 24000 --threads 2 --flip-scan --out ${out_dir}/bgen_filt_chr${i}_flip
#wait 
python -m peddy -p 4 --plot --prefix ${out_dir}/peddy_chr${i} ${out_dir}/snp_filt_chr${i}.vcf.gz ${out_dir}/snp_filt_chr${i}.ped