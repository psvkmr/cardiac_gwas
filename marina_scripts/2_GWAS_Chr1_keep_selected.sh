#!/bin/bash

#SBATCH -p shared  #batch partition
#SBATCH -J Plink_UKB    #job name
#SBATCH --output=%j.out  #output file
#SBATCH --mem 120G   #allocated memory in MB if not specified
#SBATCH -n 1	   #2 node, you can increase it



#Load Plink  module
module load apps/plink2/2.0.0a2

#This is the initial data input
plink2 --pfile plink2 --keep-fam IDs_final.txt --make-pgen --out GWAS_QC1_Chr1
