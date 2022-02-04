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


Rscript ${saige}/step2_SPAtests.R \
        --vcfFile=${out_dir}/all_chr_snpfilt_ds_only.vcf.gz \
        --vcfFileIndex=${out_dir}/all_chr_snpfilt_ds_only.vcf.gz.csi \
        --sampleFile=${out_dir}/all_chr_snpfilt.sample \
        --vcfField=DS \
        --chrom=1 \
        --minMAC=1 \
        --GMMATmodelFile=${out_dir}/SAIGE_step1.rda \
        --varianceRatioFile=${out_dir}/SAIGE_step1.varianceRatio.txt \
        --SAIGEOutputFile=${out_dir}/SAIGE_step2_output.txt \
        --numLinesOutput=2 \
        --LOCO=FALSE \
        &> ${out_dir}/SAIGE_step2.out
