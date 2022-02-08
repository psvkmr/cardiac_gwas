#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=saige_prep
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/array/saige_prep_%A_%a.out
#SBATCH --array=1-2

# load plink2
module load apps/plink2/2.0.0a2

# set variable paths for files
base_dir=/scratch/users/k2142172
out_dir=${base_dir}/outputs/cardiac_gwas/gwas_run
env=${base_dir}/packages/anaconda3/envs/RSAIGE

# set variable paths for tools
saige=${base_dir}/packages/SAIGE
bgzip=${base_dir}/packages/anaconda3/envs/peddy/bin/bgzip
bcftools=${base_dir}/packages/bcftools-1.14/bcftools

# activate SAIGE environment
conda activate $env

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/snpfilt_chr${i}*


# assign first ID column of final fam file as new .sample file for SAIGE input
#awk '{print $1}' ${out_dir}/snpfilt_chr${i}.fam > ${out_dir}/snpfilt_chr${i}.sample
#wait

# submit R script which merges res_distensibility phenotype with PCA PC values as one phenotype file
#Rscript ${base_dir}/scripts/guys_projects/cardiac_gwas/make_pheno.R ${i}
#wait

# compress the VCF file version of final data files
$bgzip -c ${out_dir}/snpfilt_chr${i}.vcf > ${out_dir}/snpfilt_chr${i}.vcf.gz
wait

# index created compressed VCF file
$bcftools index ${out_dir}/snpfilt_chr${i}.vcf.gz
wait

# SAIGE requires FORMAT to contain dosage DS info only, no genotypes GT
# so remove GT info and re-index
$bcftools annotate -x FORMAT/GT ${out_dir}/snpfilt_chr${i}.vcf.gz -O z -o ${out_dir}/snpfilt_ds_chr${i}.vcf.gz
$bcftools index ${out_dir}/snpfilt_ds_chr${i}.vcf.gz

