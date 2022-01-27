#!/bin/bash

#SBATCH -p shared  #batch partition
#SBATCH -J Plink_UKB    #job name
#SBATCH --output=%j.out  #output file
#SBATCH --mem 50G   #allocated memory in MB if not specified
#SBATCH -n 1       #2 node, you can increase it



#Load Plink  module
module load apps/plink2/2.0.0a2

#This is the final QC 
plink2 --pfile GWAS_QC2_Chr1 --maf 0.01 --geno 0.10 --hwe 0.0000000001 --make-pgen --out GWAS_QC3_Chr1
