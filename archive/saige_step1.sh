#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=saige_1
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/saige_step1.out

# load plink2
module load apps/plink2/2.0.0a2

# set variable paths for files
base_dir=/scratch/users/k2142172
out_dir=${base_dir}/outputs/cardiac_gwas/gwas_run
env=${base_dir}/packages/anaconda3/envs/RSAIGE

# set variable paths for tools
saige=${base_dir}/packages/SAIGE

# activate SAIGE environment
conda activate $env


Rscript ${saige}/step1_fitNULLGLMM.R \
        --plinkFile=${out_dir}/all_chr_snpfilt \
        --phenoFile=${out_dir}/all_chr_pca.pheno \
        --phenoCol=res_distensibility \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --outputPrefix=${out_dir}/SAIGE_step1 \
        --nThreads=2 \
        &> ${out_dir}/SAIGE_step1.out
