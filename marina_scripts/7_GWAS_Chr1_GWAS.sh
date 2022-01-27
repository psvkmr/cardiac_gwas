#!/bin/bash

#SBATCH -p shared  #batch partition
#SBATCH -J Plink_UKB    #job name
#SBATCH --output=%j.out  #output file
#SBATCH --mem 50G   #allocated memory in MB if not specified
#SBATCH -n 1       #2 node, you can increase it



#Load Plink  module
module load apps/plink2/2.0.0a2

#This is the final QC 
plink2 \
--pfile GWAS_QC3_Chr1 \
--pheno /scratch/users/stwb3495/Data/UKBiobank/res_distensibility3.txt \
--covar plink2.eigenvec \
--pheno-name res_distensibility \
--covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
--glm hide-covar \
--ci 0.95 \
--out GWAS_Results_Chr1_2
