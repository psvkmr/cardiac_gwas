#!/bin/bash -l

#SBATCH --partition=brc,shared
#SBATCH --job-name=dprep
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/dis/saige_prep_%A_%a.out
#SBATCH --array=[1-22]%6

# load plink2
module load apps/plink2/2.0.0a2

# set variable paths for files
base_dir=/scratch/users/k2142172
out_dir=${base_dir}/outputs/cardiac_gwas/dis_run
env=${base_dir}/packages/anaconda3/envs/RSAIGE

# set variable paths for tools
saige=${base_dir}/packages/SAIGE
bgzip=${base_dir}/packages/anaconda3/envs/peddy/bin/bgzip
bcftools=${base_dir}/packages/bcftools-1.14/bcftools

# activate SAIGE environment
conda activate $env

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/sample_filt_chr${i}*

# assign first ID column of final fam file as new .sample file for SAIGE input
if [ $i == 1 ]
then
awk '{print $1}' ${out_dir}/sample_filt_chr${i}.fam > ${out_dir}/sample_filt.sample
Rscript ${base_dir}/scripts/guys_projects/cardiac_gwas/dis/make_pheno.R
fi
wait


# SAIGE requires FORMAT to contain dosage DS info only, no genotypes GT
# so remove GT info and re-index
$bcftools annotate -x FORMAT/GT ${out_dir}/sample_filt_chr${i}.vcf.gz -O z -o ${out_dir}/sample_filt_ds_chr${i}.vcf.gz
wait
$bcftools index ${out_dir}/sample_filt_ds_chr${i}.vcf.gz

#sbatch prep_saige.sh /scratch/users/k2142172/outputs/cardiac_gwas/sample_ids/min_aortic_area_phenotype.txt /scratch/users/k2142172/outputs/cardiac_gwas/min_run/pca.eigenvec
