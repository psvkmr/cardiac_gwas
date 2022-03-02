run order:

filter_bgen.sh
filter_snps.sh
qc_samples.sh
filter_samples.sh <qcd_samples_to_remove.txt>
prune_snps.sh reformats.sh
{merge_chr.sh}
pca.sh {pca_all.sh}
plink_gwas.sh prep_saige.sh
step1_saige.sh
step2_saige.sh

