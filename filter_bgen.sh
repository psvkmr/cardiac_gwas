#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=min_bgen
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/min/bgen_%A_%a.out
#SBATCH --array=[1-22]%6


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

# links 
#https://github.com/alanmichaelpittman100/GWAS-Imputation-Protocol-/blob/master/HRC-1000G-check-bim-v4.2_AP_MODIFIED.pl 
#https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)

#Load Plink module
module load apps/plink2/2.0.0a2

sample_file=/scratch/users/k2142172/outputs/cardiac_gwas/sample_ids/ukb_sqc_eur_minIDs_IID.txt
raw_data=/scratch/datasets/ukbiobank/June2017/Imputed
marina_data=/scratch/users/stwb3495/Shared_Folder_PS
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/min_run

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${raw_data}/ukb_imp_chr${i}_v3_MAF1_INFO4.bgen


plink2 --bgen ${raw_data}/ukb_imp_chr${i}_v3_MAF1_INFO4.bgen ref-first \
  	--sample ${marina_data}/ukb22828_c${i}_b0_v3_s487253.sample \
	--memory 24000 \
	--threads 2 \
	--keep-fam $sample_file \
	--make-pgen \
	--out ${out_dir}/bgen_filt_chr${i}
