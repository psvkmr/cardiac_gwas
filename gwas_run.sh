srun -p brc --mem 5G --pty /bin/bash

#Load Plink  module
module load apps/plink2/2.0.0a2

#This is the initial data input

raw_data=/scratch/datasets/ukbiobank/June2017/Imputed
out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run

mkdir -p $out_dir

for i in seq `1 22`;
do
  plink2 --bgen ${data_path}/ukb_imp_chr${i}_v3_MAF1_INFO4.bgen ref-first  \
  --sample ${out_dir}/ukb22828_c1_b0_v3_s487253.sample
