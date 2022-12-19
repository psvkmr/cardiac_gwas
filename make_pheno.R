#args=commandArgs(trailingOnly=T)

out.dir <- '/scratch/users/k2142172/outputs/cardiac_gwas/min_run/'

# add sqc file?
#f <- read.table(args[1], header=T)
#p <- read.table(args[2], header=T, comment.char='')

f <- '/scratch/users/k2142172/outputs/cardiac_gwas/sample_ids/min_aortic_area_phenotype.txt'
p <- paste0(out.dir, '/genotype_pca.eigenvec')

f <- read.table(f, header = T)
p <- read.table(p, header = T, comment.char = '')

df <- merge(f, p, by.x = c('FID', 'IID'), by.y = c('X.FID', 'IID'))

write.table(df, paste0(out.dir, 'pca_pheno.pheno', sep = ''), row.names = F, quote = F, sep = '\t')
