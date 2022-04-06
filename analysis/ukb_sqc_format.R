
base.dir <- 'C:/Users/Prasanth/Documents/cardiac_gwas/'

# full all ukb samples fam file  
fam <- read.table(paste0(base.dir, 'ukb_data/ukb22418_c3_b0_v2_s488232.fam'))
# all ukb samples metadata
sqc <- read.table(paste0(base.dir, 'ukb_data/ukb_sqc_v2.txt'))
# distensibility phenotype
dis.pheno <- read.table(paste0(base.dir, 'sample_data/res_distensibility3.txt'), header = T)
# min aorta pheno
min.pheno <- read.table(paste0(base.dir, 'sample_data/min_aortic_area_phenotype.txt'), header = T)

# need to reformat metadata
sqc.file.names <- readr::read_table('ukb_data/ukb_sqc_v2_fields.txt', skip = 4, col_names = F)

sqc.names <- data.frame(col_id = sqc.file.names$X1, 
                        col_name = sqc.file.names$X2, 
                        col_class = sqc.file.names$X3, 
                        col_desc = apply(sqc.file.names[, 4:ncol(sqc.file.names)], 1, function(x) paste(x, collapse = ' ')))

sqc.names$col_name <- gsub('\\.', '_', sqc.names$col_name)

sqc.names.vec <- c('ID1', 'ID2', sqc.names[1:23, 2], paste0(rep('PC', 40), 1:40), sqc.names[25:27, 2])                        

names(sqc) <- sqc.names.vec
names(fam) <- c('FID', 'IID', 'PID', 'MID', 'SEX', 'PHENO')

sqc.full <- cbind(fam, sqc)

# distensibility phenotype merge
dis.pheno.subset <- dis.pheno[, -2]
sqc.dis <- merge(sqc.full, dis.pheno.subset, 'FID')

# min aorta phenotype merge
min.pheno.subset <- min.pheno[, -2]
sqc.min <- merge(sqc.full, min.pheno.subset, 'FID')

# distensibility subset of sqc and pheno
sqc.dis.eur <- sqc.dis[sqc.dis$in_white_British_ancestry_subset == 1, ]
# min aortic size subset
sqc.min.eur <- sqc.min[sqc.min$in_white_British_ancestry_subset == 1, ]

# add phenotype data 

#write files
# write.table(sqc.full, paste0(base.dir, 'sample_data/ukb_sqc_wIDs.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
# write.table(sqc.dis, paste0(base.dir, 'sample_data/ukb_sqc_dis.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
# write.table(sqc.min, paste0(base.dir, 'sample_data/ukb_sqc_min.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
# write.table(sqc.dis.eur, paste0(base.dir, 'sample_data/ukb_sqc_eur_dis.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
# write.table(sqc.min.eur, paste0(base.dir, 'sample_data/ukb_sqc_eur_min.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
# write.table(sqc.dis[2], paste0(base.dir, 'sample_data/ukb_sqc_dis_IID.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
# write.table(sqc.min[2], paste0(base.dir, 'sample_data/ukb_sqc_min_IID.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
# write.table(sqc.dis.eur[2], paste0(base.dir, 'sample_data/ukb_sqc_eur_dis_IID.txt'), row.names = F, col.names = T, sep = '\t', quote = F)
# write.table(sqc.min.eur[2], paste0(base.dir, 'sample_data/ukb_sqc_eur_min_IID.txt'), row.names = F, col.names = T, sep = '\t', quote = F)


#explore 'variable'
sqc.full[1:100, c('PHENO', 'Batch')]
identical(sqc.full$PHENO, sqc.full$Batch)
# The fam PHENO column contains 'redacted' factor which is missing in sqc
cbind(as.data.frame(table(sqc.full$Batch)),
as.data.frame(table(sqc.full$PHENO)))

# missing IDs from distensibility IDs not in fam are not redacted samples
sum(ids$IID %in% fam[fam$PHENO == 'redacted3', ]$V2)
# redacted samples have FID, IID set to -Int anyway
