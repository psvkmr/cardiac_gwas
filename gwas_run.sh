#!/bin/bash

#SBATCH --partition=brc
#SBATCH --time=04:00:00
#SBATCH --mem=30G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=gwas
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/gwas.out

# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

# links
#https://github.com/alanmichaelpittman100/GWAS-Imputation-Protocol-/blob/master/HRC-1000G-check-bim-v4.2_AP_MODIFIED.pl
#https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)


#Load Plink  module
module load apps/plink2/2.0.0a2

#This is the initial data input

#raw_data=/scratch/datasets/ukbiobank/June2017/Imputed
raw_data=/scratch/users/k2142172/outputs/cardiac_gwas/ukb_files
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run

mkdir -p $out_dir

for i in seq `1 2`;
do
  plink2 --bgen ${raw_data}/ukb_imp_chr${i}_v3_MAF1_INFO4.bgen ref-first  \
  --sample ${raw_data}/ukb22828_c${i}_b0_v3_s487253.sample
done

# activate SAIGE
# conda activate RSAIGE
#  FLAGPATH=`which python | sed 's|/bin/python$||'`
#  export LDFLAGS="-L${FLAGPATH}/lib"
#  export CPPFLAGS="-I${FLAGPATH}/include"
