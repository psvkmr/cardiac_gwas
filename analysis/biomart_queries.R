library(biomaRt)

mart <- useMart('ENSEMBL_MART_SNP')
ds <- useDataset('hsapiens_snp', mart)

saige.sig.snp.query <- 
getBM(mart = ds, 
      attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 
                     'chrom_strand', 'allele', 'validated', 'allele_1', 
                     'minor_allele', 'minor_allele_freq', 'clinical_significance', 
                     'study_description', 'associated_gene', 'phenotype_name', 
                     'phenotype_description', 'associated_variant_risk_allele', 
                     'ensembl_gene_name'), 
      filters = 'snp_filter', 
      values = list(saige.sig$SNPID))

saige.sig.snp.info <- merge(saige.sig, saige.sig.snp.query, by.x = 'SNPID', by.y = 'refsnp_id')
View(saige.sig.snp.info[, c('SNPID', 'Allele1', 'Allele2', 'AC_Allele2', 'AF_Allele2', 
                            'allele', 'allele_1', 'minor_allele', 'minor_allele_freq')])


###########################################################################
#VEP

# api doesn't work
# library(httr)
# library(jsonlite)
# library(xml2)
# 
# server <- "https://rest.ensembl.org"
# ext <- "/variation/homo_sapiens"
# ids <- saige.sig$SNPID[1:50]
# r <- POST(paste(server, ext, sep = ""), content_type("application/json"), accept("application/json"), body = '{ "ids" : ["rs56116432", "COSM476" ], "phenotypes" : 1, "pops" : 1, "population_genotypes" : 1 }')
# r <- POST(paste(server, ext, sep = ""), body = toJSON(list("ids"=ids, "phenotypes"=1, "pops"=1, "population_genotypes"=1)))
# # use this if you get a simple nested list back, otherwise inspect its structure
# # head(data.frame(t(sapply(content(r),c))))
# toJSON(content(r))

write.table(saige.sig$SNPID, 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/saige/queries/rsids_for_vep.txt', 
            col.names = F, row.names = F, quote = F)


#############################################################################
