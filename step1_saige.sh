#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=saige_1
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/saige_step1_%A_%a.out
#SBATCH --array=[1-22]%6

# load plink2
module load apps/plink2/2.0.0a2

# set variable paths for files
base_dir=/scratch/users/k2142172
out_dir=${base_dir}/outputs/cardiac_gwas/gwas_run
env=${base_dir}/packages/anaconda3/envs/RSAIGE

# set variable paths for tools
saige=${base_dir}/packages/SAIGE

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/pruned_chr${i}*

# activate SAIGE environment
conda activate $env


Rscript ${saige}/step1_fitNULLGLMM.R \
        --plinkFile="${out_dir}/pruned_chr${i}" \
        --phenoFile=${out_dir}/chr${i}.pheno \
        --phenoCol=res_distensibility \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative \
        --invNormalize=TRUE \
        --outputPrefix=${out_dir}/SAIGE_step1_chr${i} \
        --nThreads=2 \
        --IsOverwriteVarianceRatioFile=TRUE \
        &> ${out_dir}/SAIGE_step1_chr${i}.out

