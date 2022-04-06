conda activate ldsc

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/dis_run
ldsc_scripts=/scratch/users/k2142172/packages/ldsc
ldsc_resources=/scratch/users/k2142172/resources/ldsc_resources
#######################################################################################
# for saige sumstats

#cmunge with
python ${ldsc_scripts}/munge_sumstats.py \
  --sumstats ${out_dir}/saige_summary_statistics.txt \
  --N-col N \
  --snp SNPID \
  --a1 Allele2 \
  --a2 Allele1 \
  --signed-sumstats BETA,0 \
  --p p.value \
  --info imputationInfo \
  --merge-alleles ${ldsc_resources}/eur_w_ld_chr/w_hm3.snplist \
  --chunksize 500000 \
  --out ${out_dir}/ldsc_saige_summary_statistics_munge

python ${ldsc_scripts}/ldsc.py \
--h2 ${out_dir}/ldsc_saige_summary_statistics_munge.sumstats.gz \
--ref-ld-chr ${ldsc_resources}/eur_w_ld_chr/ \
--w-ld-chr ${ldsc_resources}/eur_w_ld_chr/ \
--out ${out_dir}/ldsc_saige_summary_statistics_h2


###############################################################################################
# plink sum stats

# munge with:
# remove 'A1' as conflicts with 'ALT'
# cut -d '  ' -f 6 --complement plink_summary_statistics.txt > plink_summary_statistics_ldsc.txt
python ${ldsc_scripts}/munge_sumstats.py \
  --sumstats ${out_dir}/plink_summary_statistics_ldsc.txt \
  --N-col OBS_CT \
  --snp ID \
  --a1 ALT \
  --a2 REF \
  --signed-sumstats BETA,0 \
  --p P \
  --merge-alleles ${ldsc_resources}/eur_w_ld_chr/w_hm3.snplist \
  --chunksize 500000 \
  --out ${out_dir}/ldsc_plink_summary_statistics_munge

# analyse with
python ${ldsc_scripts}/ldsc.py \
--h2 ${out_dir}/ldsc_plink_summary_statistics_munge.sumstats.gz \
--ref-ld-chr ${ldsc_resources}/eur_w_ld_chr/ \
--w-ld-chr ${ldsc_resources}/eur_w_ld_chr/ \
--out ${out_dir}/ldsc_plink_summary_statistics_h2

###################################################################################################
# plink genetic correlation

/scratch/users/k2142172/packages/ldsc/ldsc.py \
--rg /scratch/users/k2142172/outputs/cardiac_gwas/min_run/plink_summary_statistics_ldsc.txt.sumstats.gz,/scratch/users/k2142172/outputs/cardiac_gwas/dis_run/plink_summary_statistics_ldsc.txt.sumstats.gz \
--ref-ld-chr /scratch/users/k2142172/resources/ldsc_resources/eur_w_ld_chr/ \
--w-ld-chr /scratch/users/k2142172/resources/ldsc_resources/eur_w_ld_chr/ \
--out /scratch/users/k2142172/outputs/cardiac_gwas/cardiac_ldsc/min_v_dis_plink_sumstats_ldsc


#########################################################################################################
# partitioned h2

python ${ldsc_scripts}/ldsc.py \
	--h2 ${out_dir}/ldsc_saige_summary_statistics_munge.sumstats.gz \
	--ref-ld-chr ${ldsc_resources}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. \
	--w-ld-chr ${ldsc_resources}/weights_hm3_no_hla/weights. \
	--overlap-annot \
	--frqfile-chr ${ldsc_resources}/1000G_Phase3_frq/1000G.EUR.QC. \
	--out ${out_dir}/ldsc_saige_summary_statistics_partitioned_h2


python ${ldsc_scripts}/ldsc.py \
	--h2 ${out_dir}/ldsc_saige_summary_statistics_munge.sumstats.gz \
	--ref-ld-chr ${ldsc_resources}/1000G_Phase3_cell_type_groups/cell_type_group.2.,${ldsc_resources}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. \
	--w-ld-chr ${ldsc_resources}/weights_hm3_no_hla/weights. \
	--overlap-annot \
	--frqfile-chr ${ldsc_resources}/1000G_Phase3_frq/1000G.EUR.QC. \
	--out ${out_dir}/ldsc_saige_summary_statistics_partitioned_cell_type_h2 \
  --print-coefficients


#########################################################################################################
# cell type expression partitioned h2

# working only when running from the ldsc resources dir for some reason 
python ${ldsc_scripts}/ldsc.py \
--h2-cts ${out_dir}/ldsc_saige_summary_statistics_munge.sumstats.gz \
--ref-ld-chr ${ldsc_resources}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. \
--out ${out_dir}/ldsc_saige_summary_statistics_tissue_expression_h2 \
--ref-ld-chr-cts ${ldsc_resources}/Multi_tissue_gene_expr.ldcts \
--w-ld-chr ${ldsc_resources}/weights_hm3_no_hla/weights.
