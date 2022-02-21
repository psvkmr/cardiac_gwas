# sample QC

setwd('C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas')

het.files <- lapply(X = list.files(pattern = '*.het',
                                   path = 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/', 
                                   full.names = T), 
                    read.table, 
                    comment.char = '', 
                    header = T)

het.names <- lapply(X = list.files(pattern = '*.het',
                                   path = 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/', 
                                   full.names = T),
                    function(x) gsub('^.*snp_filt_', '', x))

names(het.files) <- het.names                    

hetChecker <- function(x){
#  x$het_rate <- ((x$`O.HOM.` - x$`E.HOM.`) / x$`O.HOM`)
#  y <- subset(x, (x$het_rate < (mean(x$het_rate) - 3*sd(x$het_rate)) | x$het_rate > (mean(x$het_rate) + 3*sd(x$het_rate))))
  y <- subset(x, (F < -0.2 | F > 0.2))
  return(y$IID)
}
het.fail <- lapply(het.files, hetChecker)
het.fail.tbl <- sort(table(unlist(het.fail)), decreasing = T)

# not removing any samples for now, possible threshold for F +- 0.2
#lapply(het.files, function(x) x[x$IID == '1543065', 'F'])

# large list of all snp alt counts 
# acount.files <- lapply(X = list.files(pattern = '*.acount', 
#                                       path = 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/', 
#                                       full.names = T), 
#                        read.table, 
#                        comment.char = '', 
#                        header = T)

smiss.files <- lapply(X = list.files(pattern = '*.smiss', 
                                     path = 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/', 
                                     full.names = T), 
                      read.table, 
                      comment.char = '', 
                      header = T)

smiss.names <- lapply(X = list.files(pattern = '*.smiss',
                                   path = 'C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/', 
                                   full.names = T),
                    function(x) gsub('^.*snp_filt_', '', x))

names(smiss.files ) <- smiss.names
smiss.full <- Reduce(function(x, y) merge(x, y, by = 'IID', all = T), smiss.files)
sample.miss.rate <- rowMeans(smiss.full[, grepl('F_MISS', names(smiss.full))])
sample.miss.fail <- smiss.full[which(sample.miss.rate > 0.1), 'IID']
