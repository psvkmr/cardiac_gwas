# sample QC

#gwas.dir <- 'C:/Users/Prasanth/Documents/cardiac_gwas/min_aorta_gwas'
gwas.dir <- 'C:/Users/Prasanth/Documents/cardiac_gwas/dis_run'

# plink-calculated inbreeding coefficient and expected vs observed heterozygosity rates
het.files <- lapply(X = list.files(pattern = '*.het$',
                                   path = paste0(gwas.dir, '/processing/het'),
                                   full.names = T),
                    read.table, 
                    comment.char = '', 
                    header = T)

het.names <- lapply(X = list.files(pattern = '*.het$',
                                   path = paste0(gwas.dir, '/processing/het'),
                                   full.names = T),
                    function(x) gsub('^.*snp_filt_', '', x))

names(het.files) <- het.names                    

# hetChecker <- function(x){
# #  x$het_rate <- ((x$`O.HOM.` - x$`E.HOM.`) / x$`O.HOM`)
# #  y <- subset(x, (x$het_rate < (mean(x$het_rate) - 3*sd(x$het_rate)) | x$het_rate > (mean(x$het_rate) + 3*sd(x$het_rate))))
#   y <- subset(x, (F < -0.2 | F > 0.2))
#   return(y$IID)
# }
# het.fail <- lapply(het.files, hetChecker)
# het.fail.tbl <- sort(table(unlist(het.fail)), decreasing = T)

# get matrices of observed and expected hom rates per chromosome
o.hom.sum <- c()
e.hom.sum <- c()
for (i in 1:length(het.files)){
  o.hom <- het.files[[i]]$`O.HOM.`
  e.hom <- het.files[[i]]$`E.HOM.`
  o.hom.sum <- cbind(o.hom.sum, o.hom)
  e.hom.sum <- cbind(e.hom.sum, e.hom)
}

# bind matrices together, rename by chromosome
het.rate <- cbind(o.hom.sum, e.hom.sum)
colnames(het.rate)[1:22] <- paste0('O_HOM_', het.names)
colnames(het.rate)[23:44] <- paste0('E_HOM_', het.names)
#map.het.rate <- t(apply(het.rate, 1, function(z) mapply(function(x, y) (x -y)/x, z[1:22], z[23:44])))
#het.rate <- cbind(het.rate, map.het.rate)

# get total het rate for observed and expected, and find which indices are samples
# that give het values more than 3 sd away from mean, either high or low
o.hom.full <- rowSums(het.rate[, 1:22])
e.hom.full <- rowSums(het.rate[, 23:44])
het.rate.full <- (o.hom.full - e.hom.full) / o.hom.full
het.rate.high <- which(het.rate.full > (mean(het.rate.full) + (3*sd(het.rate.full))))
het.rate.low <- which(het.rate.full < (mean(het.rate.full) - (3*sd(het.rate.full))))
# write FID and IID df to remove from plink dataset
het.outliers.df <- het.files$chr1_het.het[c(het.rate.high, het.rate.low), c('X.FID', 'IID')]

#########################################################################################

# large list of all snp alt counts 
# acount.files <- lapply(X = list.files(pattern = '*.acount', 
#                                       path = 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/', 
#                                       full.names = T), 
#                        read.table, 
#                        comment.char = '', 
#                        header = T)

#######################################################################################

# n. of missing genotypes per sample per chr
smiss.files <- lapply(X = list.files(pattern = '*.smiss$',
                                     path = paste0(gwas.dir, '/processing/miss'),
                                     full.names = T),
                      read.table, 
                      comment.char = '', 
                      header = T)

smiss.names <- lapply(X = list.files(pattern = '*.smiss$',
                                     path = paste0(gwas.dir, '/processing/miss'),
                                     full.names = T),
                    function(x) gsub('^.*snp_filt_', '', x))

names(smiss.files ) <- smiss.names

# merge all chr files by IID
smiss.full <- Reduce(function(x, y) merge(x, y, by = 'IID', all = T), smiss.files)
# subset smiss by only missing counts
sample.miss.sum <- smiss.full[, grepl('MISSING_CT', names(smiss.full))]
# get one number of all snps in analysis in any one sample
sample.all.snps.count <- rowSums(smiss.full[1, grepl('OBS_CT', names(smiss.full))])
# calculate the number of genotypes missing vs all snps number, per sample
sample.miss.rate <- rowSums(sample.miss.sum)/sample.all.snps.count
sample.miss.fail <- smiss.full[which(sample.miss.rate > 0.05), ]
# write df with FID and IID of samples with > 5% genotypes missing
sample.miss.df <- sample.miss.fail[, c('X.FID.x', 'IID')]
names(sample.miss.df) <- c('X.FID', 'IID')

# merged df of either het or missingness fails
het.or.smiss.outliers.df <- unique.data.frame(rbind(het.outliers.df, sample.miss.df))

# now need to remove withdrawn participants
withdrawn.consent <- read.csv('C:/Users/Prasanth/Documents/cardiac_gwas/sample_data/participants_withdrawn.csv', header = F)
withdrawn.consent.samples <- het.files$chr1_het.het[het.files$chr1_het.het$IID %in% withdrawn.consent$V1, 'IID']
withdrawn.consent.samples.df <- data.frame('X.FID' = withdrawn.consent.samples, 
                                           'IID' = withdrawn.consent.samples)

het.smiss.withdrawn.df <- unique(data.frame(rbind(het.or.smiss.outliers.df, withdrawn.consent.samples.df)))

# final outliers write out file
# write.table(het.outliers.df, paste0(gwas.dir, '/processing/het_outlier_samples.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
# write.table(sample.miss.df, paste0(gwas.dir, '/processing/excess_missing_samples.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
# write.table(withdrawn.consent.samples.df, paste0(gwas.dir, '/processing/participant_withdrawn_samples.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
# write.table(het.smiss.withdrawn.df, paste0(gwas.dir, '/processing/het_smiss_withdrawn.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
