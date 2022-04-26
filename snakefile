# GWAS snakefile
# outputs required
# bgen to pfiles filtered for cardiac samples
# cardiac samples filtered out loq qual snps
# samples filtered for non european, het outlier, withrawn participants
# filter genotyped snps
# prune sample snps for pca
# reformat files for saige
# merge chromosome files for pca
# pca from genotypes
# pca from merged
# run plink gwas with either merged pca or geno pca
# run saige prep
# saige step 1
# saige step 2

#snakemake ${out_dir}/bgen_filt_chr{21..22}.pgen --cluster 'sbatch -p brc,shared -J filter_bgen --time 04:00:00 --mem 30G --ntasks 2 --cpus-per-task 8 -o filter_bgen.out --verbose' -j 16

base_dir='/scratch/users/k2142172'
out_dir=f'{base_dir}/outputs/cardiac_gwas/gwas_run'

plink2=f'{base_dir}/packages/plink2'
bgzip=f'{base_dir}/packages/anaconda3/envs/peddy/bin/bgzip'
bcftools=f'{base_dir}/packages/bcftools-1.14/bcftools'
saige_step1=f'{base_dir}/packages/SAIGE/step1_fitNULLGLMM.R'
saige_step2=f'{base_dir}/packages/SAIGE/step2_SPAtests.R'
ldsc=f'{base_dir}/packages/ldsc'

pheno_file='/scratch/users/stwb3495/Shared_Folder_PS/res_distensibility3.txt'
pheno_name='res_distensibility'
bgen_keep_samples=f'{base_dir}/outputs/cardiac_gwas/sample_ids/ukb_sqc_eur_dis_IID.txt'
remove_samples=f'{out_dir}/het_smiss_withdrawn.txt'

chr=range(21,23,1)

rule all:
    input:
#        expand('{out_dir}/bgen_filt_chr{chr}.pgen', out_dir=out_dir, chr=chr, pheno_name=pheno_name)
        expand('{out_dir}/gwas_results_chr{chr}.{pheno_name}.glm.linear', out_dir=out_dir, chr=chr, pheno_name=pheno_name)
    resources:
        t=2

rule filter_bgen:
    input:
        bgen='/scratch/datasets/ukbiobank/June2017/Imputed/ukb_imp_chr{chr}_v3_MAF1_INFO4.bgen',
        sample='/scratch/users/stwb3495/Shared_Folder_PS/ukb22828_c{chr}_b0_v3_s487253.sample',
        ids={bgen_keep_samples}
    output:
        out_pgen='{out_dir}/bgen_filt_chr{chr}.pgen',
    params:
        out_pref='{out_dir}/bgen_filt_chr{chr}',
        plink2=plink2
    resources:
        p='brc,shared',
        J='filter_bgen',
        t='04:00:00',
        mem='30G',
        n=2,
        c=8,
        o='{bsae_dir}/tests/snake/filter_bgen.out'
    log:
        '{out_dir}/filter_bgen_{chr}.log'
    shell:
        '''
        {params.plink2} --bgen {input.bgen} ref-first --sample {input.sample} --memory 30000 \
        --keep-fam {input.ids} --make-pgen --out {params.out_pref}
        '''

rule filter_snps:
    input:
        pgen=rules.filter_bgen.output,
        snps='{base_dir}/resources/rsIDs_In_HRC_Or_Genotyped.txt'
    output:
        pgen='{out_dir}/snp_filt_chr{chr}.pgen',
        frq='{out_dir}/snp_filt_chr{chr}.acount'
    params:
        in_pref='{out_dir}/bgen_filt_chr{chr}',
        out_pref='{out_dir}/snp_filt_chr{chr}',
        plink2=plink2
    resources:
        p='brc,shared',
        J='filter_snps',
        t='02:00:00',
        mem='24G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/filter_snps.out'
    log:
        '{out_dir}/filter_snps_{chr}.log'
    shell:
        '''
        {params.plink2} --pfile {params.in_pref} --memory 16000 --maf 0.01 --geno 0.05 \
        --hwe 1e-6 --snps-only --extract {input.snps} --make-pgen --out {params.out_pref}
        {params.plink2} --pfile {params.out_pref} --memory 16000 --freq counts --out {params.out_pref}
        '''

rule qc_samples:
    input:
        pgen=rules.filter_snps.output.pgen
    output:
        het='{out_dir}/snp_filt_chr{chr}.het',
        smiss='{out_dir}/snp_filt_chr{chr}.smiss',
        vmiss='{out_dir}/snp_filt_chr{chr}.vmiss'
    params:
        pref='{out_dir}/snp_filt_chr{chr}',
        plink2=plink2
    resources:
        p='brc,shared',
        J='qc_samples',
        t='02:00:00',
        mem='24G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/qc_samples.out'
    log:
        '{out_dir}/qc_samples_{chr}.log'
    shell:
        '''
        {params.plink2} --pfile {params.pref} --memory 24000 --het --out {params.pref}
        wait
        {params.plink2} --pfile {params.pref} --memory 24000 --missing --out {params.pref}
        '''

rule filter_samples:
    input:
        pgen=rules.filter_snps.output.pgen,
        miss={remove_samples}
    output:
        out_pgen='{out_dir}/sample_filt_chr{chr}.pgen',
        out_psam='{out_dir}/sample_filt_chr{chr}.psam'
    params:
        in_pref='{out_dir}/snp_filt_chr{chr}',
        out_pref='{out_dir}/sample_filt_chr{chr}',
        plink2=plink2
    resources:
        p='brc,shared',
        J='filter_samples',
        t='02:00:00',
        mem='24G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/filter_samples.out'
    log:
        '{out_dir}/filter_samples_{chr}.log'
    shell:
        '''
        {params.plink2} --pfile {params.in_pref} --memory 10000 --remove {input.miss} \
        --make-pgen --out {params.out_pref}
        '''

rule linkage_impute:
    input:
        pgen=rules.filter_samples.output.out_pgen
    output:
        prune='{out_dir}/sample_filt_chr{chr}.prune.in'
    params:
        pref='{out_dir}/sample_filt_chr{chr}',
        plink2=plink2
    resources:
        p='brc,shared',
        J='linkage_impute',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/linkage_impute.out'
    log:
        '{out_dir}/linkage_impute_{chr}.log'
    shell:
        '''
        {params.plink2} --pfile {params.pref} --memory 5000 \
        --rm-dup exclude-mismatch --indep-pairwise 500 50 0.2 --out {params.pref}
        '''

rule prune_impute_pfile:
    input:
        pgen=rules.filter_samples.output.out_pgen,
        prune='{out_dir}/sample_filt_chr{chr}.prune.in',
        ldr='{base_dir}/resources/highLDRegionsIDs_GRCh37.txt'
    output:
        out_pgen='{out_dir}/pruned_chr{chr}.pgen',
    params:
        in_pref='{out_dir}/sample_filt_chr{chr}',
        out_pref='{out_dir}/pruned_chr{chr}',
        plink2=plink2
    resources:
        p='brc,shared',
        J='prune_impute_pfile',
        t='02:00:00',
        mem='24G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/prune_impute_pfile.out'
    log:
        '{out_dir}/prune_impute_pfile_{chr}.log'
    shell:
        '''
        {params.plink2} --pfile {params.in_pref} \
        --memory 10000 --exclude bed0 {input.ldr} \
        --extract {input.prune} --make-pgen --out {params.out_pref}
        '''

rule prune_impute_bfile:
    input:
        pgen=rules.prune_impute_pfile.output.out_pgen
    output:
        bed='{out_dir}/pruned_chr{chr}.bed'
    params:
        pref='{out_dir}/pruned_chr{chr}',
        plink2=plink2
    resources:
        p='brc,shared',
        J='prune_impute_bfile',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/prune_impute_bfile.out'
    log:
        '{out_dir}/prune_impute_bfile_{chr}.log'
    shell:
        '''
        {params.plink2} --pfile {params.pref} \
        --memory 10000 --make-bed --out {params.pref}
        '''

rule merge_chr:
    input:
        pgen='{out_dir}/pruned_chr21.bed'
    output:
        merge='{out_dir}/pruned_to_merge.txt',
        pgen='{out_dir}/all_pruned_merged.pgen'
    params:
        out_pref='{out_dir}/all_pruned_merged',
        pruned='{out_dir}/pruned_chr*.pgen',
        plink2=plink2
    resources:
        p='brc,shared',
        J='merge_chr',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/merge_chr.out'
    log:
        '{out_dir}/merge_chr.log'
    shell:
        '''
        ls {params.pruned} | cut -d. -f 1 > {output.merge}
        wait
        {params.plink2} --pmerge-list {output.merge} pfile --memory 10000 \
        --merge-max-allele-ct 2 --out {params.out_pref}
        wait
        {params.plink2} --pfile {params.out_pref} --memory 10000 --make-bed \
        --out {params.out_pref}
        '''

rule match_samples_HRC_snps:
    input:
        samples='{out_dir}/sample_filt_chr22.psam',
        ukbb='{base_dir}/resources/ukb_genotyped/ukb_binary_v2.bed',
        snps='{base_dir}/resources/rsIDs_In_HRC_Or_Genotyped.txt'
    output:
        out_pgen='{out_dir}/genotype_filt.pgen'
    params:
        in_pref='{base_dir}/resources/ukb_genotyped/ukb_binary_v2',
        out_pref='{out_dir}/genotype_filt',
        plink2=plink2
    resources:
        p='brc,shared',
        J='match_samples_HRC_snps',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/match_samples_HRC_snps.out'
    log:
        '{out_dir}/match_samples_HRC_snps.log'
    shell:
        '''
        {params.plink2} --bfile {params.in_pref} \
        --memory 16000 --keep {input.samples} --extract {input.snps} --make-pgen \
        --out {params.out_pref}
        '''

rule linkage_geno:
    input:
        pgen=rules.match_samples_HRC_snps.output.out_pgen
    output:
        out_pgen='{out_dir}/genotype_filt.prune.in'
    params:
        pref='{out_dir}/genotype_filt',
        plink2=plink2
    resources:
        p='brc,shared',
        J='linkage_geno',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/linkage_geno.out'
    log:
        '{out_dir}/linkage_geno.log'
    shell:
        '''
        {params.plink2} --pfile {params.pref} \
        --memory 5000 --rm-dup exclude-mismatch --indep-pairwise 500 50 0.2 \
        --out {params.pref}
        '''

rule prune_geno_pfile:
    input:
        pgen=rules.linkage_geno.output.out_pgen,
        prune='{out_dir}/genotype_filt.prune.in',
        ldr='{base_dir}/resources/highLDRegionsIDs_GRCh37.txt'
    output:
        out_pgen='{out_dir}/genotype_pruned.pgen',
    params:
        in_pref='{out_dir}/genotype_filt',
        out_pref='{out_dir}/genotype_pruned',
        plink2=plink2
    resources:
        p='brc,shared',
        J='prune_geno_pfile',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/prune_geno_pfile.out'
    log:
        '{out_dir}/prune_geno_pfile.log'
    shell:
        '''
        {params.plink2} --pfile {params.in_pref} \
        --memory 5000 --exclude bed0 {input.ldr} \
        --extract {input.prune} --make-pgen --out {params.out_pref}
        '''

rule prune_geno_bfile:
    input:
        pgen=rules.prune_geno_pfile.output.out_pgen
    output:
        bed='{out_dir}/genotype_pruned.bed'
    params:
        pref='{out_dir}/genotype_pruned',
        plink2=plink2
    resources:
        p='brc,shared',
        J='prune_geno_bfile',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/prune_geno_bfile.out'
    log:
        '{out_dir}/prune_geno_bfile.log'
    shell:
        '''
        {params.plink2} --pfile {params.pref} \
        --memory 5000 --make-bed --out {params.pref}
        '''

rule pca:
    input:
        pgen=rules.prune_geno_pfile.output.out_pgen
    output:
        eigenvec='{out_dir}/genotype_pca.eigenvec'
    params:
        in_pref='{out_dir}/genotype_pruned',
        out_pref='{out_dir}/genotype_pca',
        plink2=plink2
    resources:
        p='brc,shared',
        J='pca',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/pca.out'
    log:
        '{out_dir}/pca.log'
    shell:
        '''
        {params.plink2} --pfile {params.in_pref} --memory 40000 --pca 20 --out {params.out_pref}
        '''

rule plink_gwas:
    input:
        pgen=rules.filter_samples.output.out_pgen,
        pheno={pheno_file},
        pca=rules.pca.output.eigenvec
    output:
        gwas='{out_dir}/gwas_results_chr{chr}.{pheno_name}.glm.linear'
    params:
        in_pref='{out_dir}/sample_filt_chr{chr}',
        out_pref='{out_dir}/gwas_results_chr{chr}',
        plink2=plink2,
        pheno_name={pheno_name}
    resources:
        p='brc,shared',
        J='plink_gwas',
        t='02:00:00',
        mem='16G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/plink_gwas.out'
    log:
        '{out_dir}/plink_gwas_{pheno_name}_{chr}.log'
    shell:
        '''
        {params.plink2} --pfile {params.in_pref} --memory 10000 --covar {input.pca} \
        --pheno {input.pheno} --pheno-name {params.pheno_name} \
        --covar-name PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20 \
        --glm hide-covar --ci 0.95 --out {params.out_pref}
        '''

rule plink_to_vcf:
    input:
        pgen=rules.filter_samples.output.out_pgen
    output:
        out_bed='{out_dir}/sample_filt_chr{chr}.bed',
        vcf='{out_dir}/sample_filt_chr{chr}.vcf.gz',
        vcfi='{out_dir}/sample_filt_chr{chr}.vcf.gz.gzi'
    params:
        pref='{out_dir}/sample_filt_chr{chr}',
        vcf_tmp='{out_dir}/sample_filt_chr{chr}.vcf',
        plink2=plink2,
        bgzip=bgzip
    resources:
        p='brc,shared',
        J='plink_to_vcf',
        t='04:00:00',
        mem='40G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/plink_to_vcf.out'
    log:
        '{out_dir}/plink_to_vcf_{chr}.log'
    shell:
        '''
        {params.plink2} --pfile {params.pref} --make-bed --out {params.pref}
        {params.plink2} --pfile {params.pref} --export vcf vcf-dosage=DS-force \
        --out {params.pref}
        {params.bgzip} -i {params.vcf_tmp}
        '''

rule dosage_vcf:
    input:
        vcf=rules.plink_to_vcf.output.vcf
    output:
        ds_vcf='{out_dir}/sample_filt_ds_chr{chr}.vcf.gz',
        ds_vcfi='{out_dir}/sample_filt_ds_chr{chr}.vcf.gz.csi'
    params:
        bcftools=bcftools
    resources:
        p='brc,shared',
        J='dosage_vcf',
        t='04:00:00',
        mem='32G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/dosage_vcf.out'
    log:
        '{out_dir}/dosage_vcf_{chr}.log'
    shell:
        '''
        {params.bcftools} annotate -x FORMAT/GT {input.vcf} -O z -o {output.ds_vcf}
        {params.bcftools} index {output.ds_vcf}
        '''

# change to use snakemake within R with S4 object
# fix outdir variable in Rscript, not outputting to outdir atm
rule saige_sample_pheno:
    input:
        psam='{out_dir}/sample_filt_chr22.psam',
        pheno={pheno_file},
        pca=rules.pca.output.eigenvec,
    output:
        sample='{out_dir}/sample_filt.sample',
        pheno_pca='{out_dir}/pheno_pca.pheno'
    params:
        out_dir=out_dir
    resources:
        p='brc,shared',
        J='saige_sample_pheno',
        t='00:30:00',
        mem='1G',
        o='{base_dir}/tests/snake/saige_sample_pheno.out'
    log:
        '{out_dir}/saige_sample_pheno.log'
    shell:
        '''
        awk '{{print $1}}' {input.psam} > {output.sample}
        Rscript /scratch/users/k2142172/scripts/guys_projects/cardiac_gwas/make_pheno.R \
        {input.pheno} {input.pca} {params.out_dir}
        '''

# change out to log.out
# not working
#snakemake ${out_dir}/SAIGE_step1.rda --cluster 'sbatch -p brc,shared -J ss1 --time 04:00:00 \
#--mem 12G --ntasks 2 --cpus-per-task 8 -o saige_step1.out --verbose' -j 16 --rerun-incomplete --use-conda --conda-frontend conda -n
rule saige_step1_alt:
    conda:
        '/scratch/users/k2142172/packages/anaconda3/envs/RSAIGE/RSAIGE.yaml'
    input:
        prune=rules.prune_geno_bfile.output.bed,
        pheno_pca=rules.saige_sample_pheno.output.pheno_pca
    output:
        rda='{out_dir}/SAIGE_step1.rda',
        marks='{out_dir}/SAIGE_step1_30markers.SAIGE.results.txt',
        var='{out_dir}/SAIGE_step1.varianceRatio.txt',
        out='{out_dir}/SAIGE_step1.out'
    params:
        pheno_name=pheno_name,
        in_pref='{out_dir}/genotype_pruned',
        out_pref='{out_dir}/SAIGE_step1',
        saige=saige_step1
    resources:
        p='brc,shared',
        J='saige_step1',
        t='12:00:00',
        mem='72G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/saige_step1.out'
    log:
        '{out_dir}/saige_step1.log'
    shell:
        '''
        Rscript {params.saige} --plinkFile={params.in_pref} --phenoFile={input.pheno_pca} \
        --phenoCol={params.pheno_name} --sampleIDColinphenoFile=IID --traitType=quantitative \
        --invNormalize=TRUE --outputPrefix={params.out_pref} --nThreads=4 \
        --IsOverwriteVarianceRatioFile=TRUE &> {output.out}
        '''

rule saige_step1:
    input:
        prune=rules.prune_geno_bfile.output.bed,
        pheno_pca=rules.saige_sample_pheno.output.pheno_pca
    output:
        rda='{out_dir}/SAIGE_step1.rda',
        marks='{out_dir}/SAIGE_step1_30markers.SAIGE.results.txt',
        var='{out_dir}/SAIGE_step1.varianceRatio.txt',
        out='{out_dir}/SAIGE_step1.out'
    params:
        pheno_name=pheno_name,
        out_dir=out_dir,
    resources:
        p='brc,shared',
        J='saige_step1',
        t='12:00:00',
        mem='72G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/saige_step1.out'
    log:
        '{out_dir}/saige_step1.log'
    shell:
        '''
        bash /scratch/users/k2142172/scripts/guys_projects/cardiac_gwas/step1_saige.sh {params.pheno_name} {params.out_dir}
        '''


rule saige_step2:
    input:
        vcf=rules.dosage_vcf.output.ds_vcf,
        vcfi=rules.dosage_vcf.output.ds_vcfi,
        sample=rules.saige_sample_pheno.output.sample,
        rda=rules.saige_step1.output.rda,
        var=rules.saige_step1.output.var
    output:
        res='{out_dir}/SAIGE_step2_chr{chr}.txt',
        out='{out_dir}/SAIGE_step2_chr{chr}.out'
    params:
        saige=saige_step2,
        chr=chr
    resources:
        p='brc,shared',
        J='saige_step2',
        t='04:00:00',
        mem='32G',
        n=2,
        c=8,
        o='{base_dir}/tests/snake/saige_step2.out'
    log:
        '{out_dir}/saige_step2_{chr}.log'
    shell:
        '''
        Rscript {params.saige} --vcfFile={input.vcf} --vcfFileIndex={input.vcfi} \
        --sampleFile={input.sample} --vcfField=DS --chrom={params.chr} --minMAC=1 \
        --GMMATmodelFile={input.rda} --varianceRatioFile={input.var} --numLinesOutput=2 \
        --SAIGEOutputFile={output.res} &> {output.out}
        '''

rule ldsc_munge_saige:
    input:
        ss='{out_dir}/saige_summary_statistics.txt',
        snplist='/scratch/users/k2142172/resources/ldsc_resources/eur_w_ld_hr/w_hm3.snplist'
    output:
        munged='{out_dir}/saige_summary_statistics_ldsc_munged.sumstats.gz'
    params:
        ldsc='{ldsc}/munge_sumstats.py',
        out_pref='{out_dir}/saige_summary_statistics_ldsc_munged'
    resources:
        p='brc,shared',
        J='ldsc_munge',
        t='01:00:00',
        mem='8G',
        o='{base_dir}/tests/snake/ldsc_munge_saige.out'
    log:
        '{out_dir}/ldsc_munge_saige.log'
    shell:
        '''
        {params.ldsc} --sumstats {input.ss} \
          --N-col N \
          --snp SNPID \
          --a1 Allele2 \
          --a2 Allele1 \
          --signed-sumstats BETA,0 \
          --p p.value \
          --info imputationInfo \
          --merge-alleles {input.snplist} \
          --chunksize 500000 \
          --out {params.out_pref}
        '''

rule ldsc_saige:
    input:
        ss=rules.ldsc_munge_saige.output.munged
    output:
        h2='{out_dir}/saige_summary_statistics_ldsc_h2.'
    params:
        ldsc='{ldsc}/ldsc.py',
        out_pref='{out_dir}/saige_summary_statistics_ldsc_h2',
        ref_ld='{base_dir}/resources/ldsc_resources/g1000_eur/1000G_Phase3_baselineLD_v2.2_ldscores/',
        w_ld='{base_dir}/resources/ldsc_resources/g1000_eur/1000G_Phase3_baselineLD_v2.2_ldscores/'
    resources:
        p='brc,shared',
        J='ldsc',
        t='01:00:00',
        mem='8G',
        o='{base_dir}/tests/snake/ldsc_saige.out'
    log:
        '{out_dir}/ldsc_saige.log'
    shell:
        '''
        {params.ldsc} --h2 {inputs.ss} --ref-ld-chr {params.ref_ld} --w-ld-chr {params.w_ld} --out {params.out_pref}
        '''

rule ldsc_munge_plink:
    input:
        ss='{out_dir}/plink_summary_statistics.txt',
        snplist='/scratch/users/k2142172/resources/ldsc_resources/eur_w_ld_hr/w_hm3.snplist'
    output:
        munged='{out_dir}/plink_summary_statistics_ldsc_munged.sumstats.gz'
    params:
        ldsc='{ldsc}/munge_sumstats.py',
        out_pref='{out_dir}/plink_summary_statistics_ldsc_munged'
    resources:
        p='brc,shared',
        J='ldsc_munge',
        t='01:00:00',
        mem='8G',
        o='{base_dir}/tests/snake/ldsc_munge_plink.out'
    log:
        '{out_dir}/ldsc_munge_plink.log'
    shell:
        '''
        {params.ldsc} --sumstats {input.ss} \
          --N-col OBS_CT \
          --snp ID \
          --a1 ALT \
          --a2 REF \
          --signed-sumstats BETA,0 \
          --p P \
          --merge-alleles /scratch/users/k2142172/resources/ldsc_resources/eur_w_ld_chr/w_hm3.snplist \
          --chunksize 500000 \
          --out {params.out_pref}
        '''

rule ldsc_plink:
    input:
        ss=rules.ldsc_munge_plink.output.munged
    output:
        h2='{out_dir}/plink_summary_statistics_ldsc_h2.'
    params:
        ldsc='{ldsc}/ldsc.py',
        out_pref='{out_dir}/plink_summary_statistics_ldsc_h2',
        ref_ld='{base_dir}/resources/ldsc_resources/g1000_eur/1000G_Phase3_baselineLD_v2.2_ldscores/',
        w_ld='{base_dir}/resources/ldsc_resources/g1000_eur/1000G_Phase3_baselineLD_v2.2_ldscores/'
    resources:
        p='brc,shared',
        J='ldsc',
        t='01:00:00',
        mem='8G',
        o='{base_dir}/tests/snake/ldsc_plink.out'
    log:
        '{out_dir}/ldsc_plink.log'
    shell:
        '''
        {params.ldsc} --h2 {inputs.ss} --ref-ld-chr {params.ref_ld} --w-ld-chr {params.w_ld} --out {params.out_pref}
        '''
