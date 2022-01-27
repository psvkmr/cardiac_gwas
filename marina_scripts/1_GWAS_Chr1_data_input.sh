#!/bin/bash

#SBATCH -p shared  #batch partition
#SBATCH -J Plink_UKB    #job name
#SBATCH --output=%j.out  #output file
#SBATCH --mem 100G   #allocated memory in MB if not specified
#SBATCH -n 1       #2 node, you can increase it



#Load Plink  module
module load apps/plink2/2.0.0a2

#This is the initial data input
plink2 --bgen /scratch/datasets/ukbiobank/June2017/Imputed/ukb_imp_chr1_v3_MAF1_INFO4.bgen ref-first --sample /scratch/users/stwb3495/Data/UKBiobank/ukb22828_c1_b0_v3_s487253.sample
