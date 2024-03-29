library(tidyverse)
library(PCAtools)

# plink v saige pvalues
###############################################################################
plink.merge <- plink.filt
saige.merge <- saige.filt

names(plink.merge)[c(1:5, 13, 14)] <- c('CHR', 'POS', 'SNPID', 'Allele1', 'Allele2', 'Tstat', 'p.value')

merged <- full_join(plink.merge, 
                    saige.merge, 
                    by = c('CHR', 'POS', 'SNPID'), 
                    suffix = c('.plink', '.saige'))

plot(-log10(merged$p.value.plink), -log10(merged$p.value.saige))
plot(-log10(merged$BETA.plink), -log10(merged$BETA.saige))
plot(-log10(merged$Tstat.plink), -log10(merged$Tstat.saige))

View(head(arrange(merged, p.value.plink)))


###########################################################################
# read pheno files

pheno.file.paths <- list.files(path = 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/pheno',
                               pattern = 'chr.*.pheno', full.names = T)
pheno.files <- lapply(pheno.file.paths, read.table, header = T)
pheno.names <- gsub('.pheno', '', gsub('^.*/pheno/', '', pheno.file.paths))
names(pheno.files) <- pheno.names

# plot(pheno.files$chr7$PC1, pheno.files$chr7$res_distensibility)
# plot(pheno.files$chr7$PC2, pheno.files$chr7$res_distensibility)
ggplot(pheno.files$chr7, aes(PC1, PC2)) + 
  geom_point(aes(colour = res_distensibility), alpha = 0.8) +
  scale_colour_viridis_c()

pcas <- lapply(pheno.files, function(x) as.matrix(x[, 3:22]))


###############################################################################
# extract sig snp genotypes

extract.saige.dosages <- saige.sig[, c('CHR', 'POS')]
# write.table(extract.saige.dosages, 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/extract_saige_dosages.tab',
#            quote = F, sep = '\t', col.names = F, row.names = F)
# use bcftools on rosalind to extract from dosage vcf files

# load vcf results
saige.sig.dosages.chr7 <- read.table('C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/saige/af/extracted_saige_sig_dosages_chr7.tab', 
                                     comment.char = '', header = F)
# load vcf sample names separately
vcf.samples <- read.table('C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/processing/vcf_sample_names.txt')
# create vcf headers from standard columns, plus loaded sample names
vcf.headers <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', vcf.samples$V1)
names(saige.sig.dosages.chr7) <- vcf.headers

# convert vcf to matrix of dosages, and transform with samples rows and snps columns
saige.sig.dos.7.mat <- dplyr::select(saige.sig.dosages.chr7, -c('CHROM', 'POS', 'REF', 'ALT', 'QUAL')) %>% 
  column_to_rownames('ID') %>% t()
mat.samples <- rownames(saige.sig.dos.7.mat)
mat.samples <- gsub('^.*_', '', mat.samples)
rownames(saige.sig.dos.7.mat) <- mat.samples

########################################################################################
# check genotype to dist, pc1 correlations

# load phenotype data for chromosome
pheno.7 <- pheno.files$chr7
pheno.7$FID <- as.character(pheno.7$FID)

# convert matrix back to dataframe and merge with phenotype data
saige.sig.dos.res.7 <- as.data.frame(saige.sig.dos.7.mat) %>% rownames_to_column('ID')
saige.sig.dos.res.7 <- left_join(pheno.7, saige.sig.dos.res.7, by = c('FID' = 'ID'))
# write.table(saige.sig.dos.res.7, 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/saige/af/samples_by_snps_chr7_dosage_matrix_w_phenotype.tab', 
#             row.names = F, quote = F, sep = '\t')

# compare dosage values of individual snps vs phenotype value across samples
# differences should follow the signed beta value in sumstats
# example of damaging snp
saige.sig.dos.res.7 %>% group_by(as.factor(rs2071307)) %>% summarise(mean(res_distensibility))
# example of protective snp 
saige.sig.dos.res.7 %>% group_by(as.factor(rs10224499)) %>% summarise(mean(res_distensibility))

######################################################################
# pca outliers

pca.files <- lapply(list.files(path = 'C:/Users/Prasanth/Documents/cardiac_gwas/gwas_results/pca',
                          pattern = 'pca_chr.*.eigenvec', full.names = T),
               read.table, header = F)
pca.filenames <- list.files(path = 'C:/Users/Prasanth/Documents/cardiac_gwas/gwas_results/pca',
                            pattern = 'pca_chr.*.eigenvec', full.names = T)
pca.filenames <- gsub('.eigenvec', '', gsub('^.*/pca/pca_', '', pca.filenames))
names(pca.files) <- pca.filenames

detectOutliers <- function(x) {
  mean.x <- mean(x)
  sd.x <- sd(x)
  outliers <- c()
  for (i in 1:length(x)){
    if (x[i] > mean.x + 3*sd.x) {
      outliers <- c(outliers, x[i])
    } else if (x[i] < mean.x - 3*sd.x) {
      outliers <- c(outliers, x[i])
    }
  }
  outliers.rows <- match(outliers, x)
  return(outliers.rows)
}

pca.filter.columns <- c(3:4)    
pc.outliers <- lapply(pca.files, function(x) x[, pca.filter.columns])
pc.outliers <- lapply(pc.outliers, function(x) lapply(x, detectOutliers))
pca.outlier.samples <- lapply(pc.outliers, function(x) unique(unlist(x)))
pca.all.outlier.samples <- unique(unlist(pca.outlier.samples))
pca.outlier.samples.df <- data.frame(fid = pca.all.outlier.samples, iid = pca.all.outlier.samples)
final.samples <- lapply(pca.files, `[[`, 1)
final.samples <- Reduce(intersect, final.samples)
final.samples <- final.samples[-pca.all.outlier.samples]
final.samples.df <- data.frame(fid = final.samples, iid = final.samples)
# write.table(final.samples.df, 
#             "C:/Users/Prasanth/Documents/cardiac_gwas/gwas_results/final_samples.txt", 
#             col.names = F, row.names = F, quote = F, sep = '\t')
final.pcas <- lapply(pca.files, function(x) x[x$V1 %in% final.samples, ])
final.pcas <- lapply(final.pcas, function(x) `names<-`(x, c('#FID', 'IID', paste0(rep('PC'), 1:10))))
# for (i in 1:length(final.pcas)){
#   write.table(final.pcas[[i]],
#               file = paste0('C:/Users/Prasanth/Documents/cardiac_gwas/strict_gwas/pcas/pca_',
#                             pca.filenames[i], '.eigenvec', sep = ''),
#               row.names = F, quote = F, sep = '\t')
# }


#####################################################################################
# mitchell
mitchell.cat <- read_delim('mitchell_gwas_pubmedid_22068335.tsv', delim = '\t')
mitchell.own <- read.csv('mitchell_own_table.csv')

mitchell.own$p.value <- gsub('\\?', '-', mitchell.own$p.value)
mitchell.own$p.value <- gsub(' x 10', 'e', mitchell.own$p.value)
mitchell.own$p.value <- as.numeric(mitchell.own$p.value)