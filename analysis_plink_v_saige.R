
# plink v saige
###############################################################################
plink.merge <- plink.filt
saige.merge <- saige.filt

names(plink.merge)[c(1:5, 13, 14)] <- c('CHR', 'POS', 'SNPID', 'Allele1', 'Allele2', 'Tstat', 'p.value')

merged <- full_join(plink.merge, 
                    saige.merge, 
                    by = c('CHR', 'POS', 'SNPID'), 
                    suffix = c('.plink', '.saige'))

plot(-log10(merged$p.value.plink), -log10(merged$p.value.saige))

View(head(arrange(merged, p.value.plink)))


###########################################################################
# read files

pheno.files <- lapply(list.files(path = 'C:/Users/Prasanth/Documents/cardiac_gwas/gwas_results/pheno', 
                                 pattern = 'chr.*.pheno', full.names = T), read.table, header = T)
#pheno.names <- lapply(pheno.files, function(x) paste0('chr', x$V1[1], sep = ''))
#names(pheno.files) <- pheno.names
lm(res_distensibility~PC1, pheno.files[[1]])
plot(pheno.files[[1]]$PC1, pheno.files[[1]]$res_distensibility)
plot(pheno.files[[1]]$PC2, pheno.files[[1]]$res_distensibility)

###############################################################################
# extract sig snp genotypes

saige.sig.file <- saige.sig[, c('CHR', 'POS')]
#write.table(saige.sig.file, 'C:/Users/Prasanth/Documents/cardiac_gwas/gwas_results/saige_sig.tab', 
#            quote = F, sep = '\t', col.names = F, row.names = F)

saige.sig.dosages.chr7 <- read.table('C:/Users/Prasanth/Documents/cardiac_gwas/gwas_results/saige_sig_dosages_chr7.tab', 
                                     comment.char = '', skip = 10, header = T)
saige.sig.dos.7.mat <- saige.sig.dosages.chr7[-2, ]
saige.sig.dos.7.mat <- dplyr::select(saige.sig.dos.7.mat, -c('X.CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')) %>% 
  `rownames<-`(NULL) %>%
  column_to_rownames('ID') %>%
  t()

rownames(saige.sig.dos.7.mat) <- gsub('^X.*_', '', rownames(saige.sig.dos.7.mat))
saige.sig.dos.7.mat <- saige.sig.dos.7.mat %>% as.data.frame() %>% rownames_to_column('ID')

########################################################################################
# check genotype to dist, pc1 correlations

pheno.7 <- pheno.files[[7]]
pheno.7$FID <- as.character(pheno.7$FID)
saige.sig.dos.res.7 <- left_join(pheno.7, saige.sig.dos.7.mat, by = c('FID' = 'ID'))
saige.sig.dos.res.7 <- na.omit(saige.sig.dos.res.7)

# cor(saige.sig.dos.res.7$PC1, saige.sig.dos.res.7$rs35361607)
# cor(saige.sig.dos.res.7$PC2, saige.sig.dos.res.7$rs35361607)
# cor(saige.sig.dos.res.7$res_distensibility, saige.sig.dos.res.7$rs35361607)
# cor(saige.sig.dos.res.7$PC1, saige.sig.dos.res.7$res_distensibility)

saige.sig.dos.res.7$dos.total <- rowSums(saige.sig.dos.res.7[, c(14:ncol(saige.sig.dos.res.7))])

cor(saige.sig.dos.res.7$PC1, saige.sig.dos.res.7$dos.total)
cor(saige.sig.dos.res.7$PC2, saige.sig.dos.res.7$dos.total)
cor(saige.sig.dos.res.7$res_distensibility, saige.sig.dos.res.7$dos.total)
cor(saige.sig.dos.res.7$PC1, saige.sig.dos.res.7$res_distensibility)
cor(saige.sig.dos.res.7$PC2, saige.sig.dos.res.7$res_distensibility)

plot(saige.sig.dos.res.7$res_distensibility, saige.sig.dos.res.7$dos.total)

######################################################################
# pca
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


###############################################################################

# mitchell
mitchell.cat <- read_delim('mitchell_gwas_pubmedid_22068335.tsv', delim = '\t')
mitchell.own <- read.csv('mitchell_own_table.csv')

mitchell.own$p.value <- gsub('\\?', '-', mitchell.own$p.value)
mitchell.own$p.value <- gsub(' x 10', 'e', mitchell.own$p.value)
mitchell.own$p.value <- as.numeric(mitchell.own$p.value)