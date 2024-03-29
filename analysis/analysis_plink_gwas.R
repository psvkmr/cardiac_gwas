library(tidyverse)

#https://github.com/GP2-TNC-WG/GP2-Bioinformatics-course/blob/master/Module_III.md

#plink.dir <- 'C:/Users/Prasanth/Documents/cardiac_gwas/min_run/plink/'
plink.dir <- 'C:/Users/Prasanth/Documents/cardiac_gwas/dis_run/plink/'

#plink.files <- lapply(list.files(pattern = 'gwas_geno.*.linear', 
plink.files <- lapply(list.files(pattern = 'gwas_genotype.*.linear', 
                                 path = paste0(plink.dir, 'outputs'), 
                                 full.names = T), 
                      read.table, comment.char = "", header = T)
plink.names <- lapply(plink.files, function(x) paste0('chr', x$`X.CHROM`[1], sep = ''))
names(plink.files) <- plink.names

plink <- Reduce(rbind, plink.files)
names(plink)[1] <- 'CHR'
write.table(plink, paste0(plink.dir, 'plink_summary_statistics.txt'),
          row.names = F, quote = F, sep = '\t')
R.utils::gzip(paste0(plink.dir, 'plink_summary_statistics.txt'), destname = paste0(plink.dir, 'plink_summary_statistics.txt.gz'), 
              overwrite = F, remove = F)

plink.sig <- filter(plink, P < 5E-8)
write.table(plink.sig, paste0(plink.dir, 'plink_sig_summary_statistics.csv'),
          row.names = F, quote = F, sep = '\t')

plink.log10 <- plink %>%
  mutate(log10Praw = -1*log10(P)) %>%
  mutate(log10P = ifelse(log10Praw > 40, 40, log10Praw)) %>%
  mutate(log10Plevel = ifelse(P < 5E-08, "possible", NA)) %>%
  mutate(log10Plevel = ifelse(P < 5E-09, "likely", log10Plevel))

# Reduction of the GWAS object
# This is for more efficient plotting
# This drops everything not useful in future AUC calcs estiamte in Nalls et al., 2019

plink.filt <- filter(plink.log10, log10P > 3.114074) %>%
  arrange(CHR, POS)

plink.plotting <- plink.log10 %>%
  group_by(CHR) %>%
  summarise(chr.len=max(POS)) %>%
  mutate(tot = cumsum(as.numeric(chr.len)) - chr.len) %>%
  dplyr::select(-chr.len) %>%
  left_join(plink.filt, ., by = 'CHR') %>%
#  left_join(plink.log10, ., by = 'CHR') %>%
  arrange(ordered(CHR), POS) %>%
  mutate(POScum=POS+tot)

plink.axis.df <- plink.plotting %>%
  group_by(CHR) %>%
  summarise(center = (max(POScum) + min(POScum)) / 2)

plink.manhattan <-
  ggplot(plink.plotting, aes(POScum, log10P)) +
  geom_point(aes(colour=as.factor(CHR)), alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("dark grey", "black"), 22)) +
  scale_x_continuous(label = plink.axis.df$CHR, breaks = plink.axis.df$center) +
  scale_y_continuous(expand = c(0, 0.05)) +
  theme_classic() +
  labs(x = "CHR", y = "-log10P") +
  theme(legend.position = "none") +
  geom_abline(intercept = -log10(5E-08), slope = 0, linetype = 2)
ggsave(paste0(plink.dir, 'plink_manhattan.png'),
       plot = plink.manhattan, dpi = 400)

plink.qqplotter <- function(df) {
  df <- df[!is.na(df$P), ]
  df <- mutate(df, exp.pvalues = (rank(df$P, ties.method="first")+.5)/(length(df$P)+1))
  ggplot(df, aes(-log10(exp.pvalues), -log10(P))) + 
    geom_point() +
    labs(x = "Expected logP", y = "Observed logP") + 
    geom_abline(slope = 1, intercept = 0) +
    #geom_label_repel(data = arrange(df, p.value)[1:5, ], box.padding = 0.5) +
    theme_classic() 
}

plink.qqplt <- plink.qqplotter(plink)
# ggsave('C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/plink_qqplot.png',
#        plot = plink.qqplt, dpi = 400)

plink.chisq <- qchisq(1-plink$P, 1)
plink.lambda <- median(plink.chisq) / qchisq(0.5,1)
write.table(plink.lambda, paste0(plink.dir, 'plink_lambda_value.txt'))

plink.alt.chisq <- qnorm(plink$P/2)
plink.alt.lambda <- median(plink.alt.chisq^2, na.rm=T)/qchisq(0.5,df=1)

plink.qqplt