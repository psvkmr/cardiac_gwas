run order:

filter_bgen.sh <sample_ids.txt>
filter_snps.sh
qc_samples.sh
filter_samples.sh <qcd_samples_to_remove.txt>
prune_snps.sh  &  reformats.sh  &  filter_genotypes.sh
merge_chr.sh  &  pca_genotypes.sh
pca_all.sh
plink_gwas.sh <phenotype file> <phenotype col>  &  prep_saige.sh <phenotype data file>  & plink_geno_gwas.sh
step1_saige.sh
step2_saige.sh

	
