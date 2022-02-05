
# plink v saige
###############################################################################
plink.merge <- plink
saige.merge <- saige

names(plink.merge)[c(1:5, 13, 14)] <- c('CHR', 'POS', 'SNPID', 'Allele1', 'Allele2', 'Tstat', 'p.value')

merged <- full_join(plink.merge, 
                    saige.merge, 
                    by = c('CHR', 'POS', 'SNPID'), 
                    suffix = c('.plink', '.saige'))

plot(-log10(merged$p.value.plink), -log10(merged$p.value.saige))

View(head(arrange(merged, p.value.plink)))

###############################################################################

# mitchell
mitchell.cat <- read_delim('mitchell_gwas_pubmedid_22068335.tsv', delim = '\t')
mitchell.own <- read.csv('mitchell_own_table.csv')

mitchell.own$p.value <- gsub('\\?', '-', mitchell.own$p.value)
mitchell.own$p.value <- gsub(' x 10', 'e', mitchell.own$p.value)
mitchell.own$p.value <- as.numeric(mitchell.own$p.value)
