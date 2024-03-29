#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=saige_2
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/saige_step2_%A_%a.out
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
ls -alht ${out_dir}/sample_filt_chr${i}*

# activate SAIGE environment
conda activate $env


Rscript ${saige}/step2_SPAtests.R \
        --vcfFile=${out_dir}/sample_filt_ds_chr${i}.vcf.gz \
        --vcfFileIndex=${out_dir}/sample_filt_ds_chr${i}.vcf.gz.csi \
        --sampleFile=${out_dir}/sample_filt_chr${i}.sample \
        --vcfField=DS \
        --chrom=${i} \
        --minMAC=1 \
        --GMMATmodelFile=${out_dir}/SAIGE_step1.rda \
        --varianceRatioFile=${out_dir}/SAIGE_step1.varianceRatio.txt \
        --SAIGEOutputFile=${out_dir}/SAIGE_step2_chr${i}.txt \
        --numLinesOutput=2 \
        --LOCO=FALSE \
        &> ${out_dir}/SAIGE_step2_chr${i}.out

