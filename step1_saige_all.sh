#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=m2saige_1
#SBATCH --time=12:00:00
#SBATCH --mem=72G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/min/saige_all_step1.out

# load plink2
module load apps/plink2/2.0.0a2

# set variable paths for files
#phenotype_col=$1
base_dir=/scratch/users/k2142172
out_dir=${base_dir}/outputs/cardiac_gwas/min_run
env=${base_dir}/packages/anaconda3/envs/RSAIGE

# set variable paths for tools
saige=${base_dir}/packages/SAIGE

ls -alht ${out_dir}/all_pruned_merged*

# activate SAIGE environment
conda activate $env


Rscript ${saige}/step1_fitNULLGLMM.R \
        --plinkFile="${out_dir}/all_pruned_merged" \
        --phenoFile=${out_dir}/pca_pheno.pheno \
        --phenoCol="residual_min_aortic_area_mm2" \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --outputPrefix=${out_dir}/SAIGE_step1_all \
        --nThreads=4 \
        --IsOverwriteVarianceRatioFile=TRUE \
        &> ${out_dir}/SAIGE_step1_all.out

