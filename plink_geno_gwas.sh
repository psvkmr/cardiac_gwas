#!/bin/bash

#SBATCH --partition=brc,shared
#SBATCH --job-name=mgplink
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --verbose
#SBATCH --output=/scratch/users/k2142172/tests/min/gwas_geno_plink_%A_%a.out
#SBATCH --array=1-22


# script exits if return value of a command is not zero
set -e
# print shell input lines as they are read for debugging
set -v
# prevents output redirection from overwriting existing files
#set -o noclobber

#Load Plink module
module load apps/plink2/2.0.0a2

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/min_run

# phenotype files and col name
# /scratch/users/k2142172/outputs/cardiac_gwas/res_distensibility3.txt res_distensibility
# /scratch/users/k2142172/outputs/cardiac_gwas/min_aortic_area_phenotype.txt residual_min_aortic_area_mm2
pheno_file=/scratch/users/k2142172/outputs/cardiac_gwas/sample_ids/ukb_sqc_eur_min.txt
pheno_name='residual_min_aortic_area_mm2'

mkdir -p $out_dir

echo $SLURM_ARRAY_TASK_ID
i=$SLURM_ARRAY_TASK_ID
ls -alht ${out_dir}/sample_filt_chr${i}*

plink2 --pfile ${out_dir}/sample_filt_chr${i} \
--memory 24000 \
--threads 2 \
--pheno $pheno_file \
--covar ${out_dir}/genotype_pca.eigenvec \
--pheno-name $pheno_name \
--covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
--glm hide-covar \
--ci 0.95 \
--out ${out_dir}/gwas_genotype_results_chr${i}
