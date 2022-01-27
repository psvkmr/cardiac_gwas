#!/bin/bash

#SBATCH -p shared  #batch partition
#SBATCH -J Plink_UKB    #job name
#SBATCH --output=%j.out  #output file
#SBATCH --mem 50G   #allocated memory in MB if not specified
#SBATCH -n 1       #2 node, you can increase it



#Load Plink  module
module load apps/plink2/2.0.0a2

#This is the initial data input
plink2 --pfile GWAS_QC1_Chr1 --extract plink2.prune.in --make-pgen --out GWAS_QC2_Chr1
