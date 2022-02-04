library(dplyr)

out_dir=/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run

f <- read.table(paste0(out_dir, 'res_distensibility3.txt', sep = ''), header=T)
p <- read.table(paste0(out_dir, 'all_chr_pca.eigenvec', sep = ''), header=F)

names(p) <- c('FID', 'IID', paste0(rep('PC'), seq(1, 10), sep = ''))
df <- left_join(p, f, by = c('FID', 'IID'))

write.table(df, paste0(out_dir, 'all_chr_pca.pheno', sep = ''), row.names = F, quote = F, sep = '\t')
