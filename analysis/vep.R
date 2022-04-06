library(tidyverse)

###############################################################################

#saige vep 
gwas.dir <- 'C:/Users/Prasanth/Documents/cardiac_gwas/dis_run/saige/'

vep <- read_table(paste0(gwas.dir, 'saige/saige_sig_summary_statistics_rsIDs.txt.out'), 
                        comment = '##')
saige.sig <- read_table(paste0(gwas.dir, 'saige_sig_summary_statistics.txt'))

# split vep 'extra' column into separate columns and rejoin to df
vep.extra <- strsplit(vep$Extra, ';') %>% lapply(strsplit, '=')
extractExtra <- function(Extra){
  df <- data.frame()
  for (i in 1:length(Extra)){
    dat <- lapply(Extra[[i]], `[[`, 2)
    names(dat) <- lapply(Extra[[i]], `[[`, 1)
    df <- bind_rows(df, dat)
  }
  return(df)
}

vep.extra.df <- extractExtra(vep.extra)
vep <- cbind(vep, vep.extra.df)

# subset vep for variables of interest
vep.sub <- vep[, c('#Uploaded_variation', 'Allele', 'Consequence', 'IMPACT', 'SYMBOL',
                   'Gene', 'BIOTYPE', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF')]
vep.sub <- unique.data.frame(vep.sub)

# merge vep with gwas 
saige.vep <- merge(saige.sig, vep, by.x = 'SNPID', by.y = '#Uploaded_variation')
# write.table(saige.vep, paste0(gwas.dir, 'vep/saige_vep_merge.tab'), row.names = F, 
#             quote = F, sep = '\t')

saige.vep.sub <- merge(saige.sig, vep.sub, by.x = 'SNPID', by.y = '#Uploaded_variation')
af.comparison <- saige.vep.sub[, c('SNPID', 'Allele1', 'Allele2', 'AF_Allele2', 
                                   'Allele', 'AF', 'AMR_AF', 'EUR_AF')]

# check concordance of allele freqs in gwas data vs vep 
af.correlation.scatter <- plot(af.comparison$AF_Allele2, af.comparison$EUR_AF)
#plot(af.comparison$AF_Allele2, af.comparison$AMR_AF)

table(saige.vep.sub$SYMBOL)
# write.table(as.data.frame(table(saige.vep.sub$SYMBOL)), paste0(gwas.dir, 'vep/variants_by_gene_count.txt'),
#             row.names = F, quote = F, sep = '\t')
