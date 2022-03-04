library(dplyr)

args=commandArgs(trailingOnly=T)

out_dir <- '/scratch/users/k2142172/outputs/cardiac_gwas/gwas_run/'

f <- read.table(args[1], header=T, fill=T)
p <- read.table(paste0(out_dir, 'pca.eigenvec', sep = '') header=F)

names(p) <- c('FID', 'IID', paste0(rep('PC'), seq(1, 20), sep = ''))
df <- left_join(p, f, by = c('FID', 'IID'))

write.table(df, paste0(out_dir, 'pheno_file.pheno', sep = ''), row.names = F, quote = F, sep = '\t')

