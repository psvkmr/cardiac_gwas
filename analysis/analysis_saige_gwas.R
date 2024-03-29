library(tidyverse)

#https://github.com/GP2-TNC-WG/GP2-Bioinformatics-course/blob/master/Module_III.md
saige.dir <- 'C:/Users/Prasanth/Documents/cardiac_gwas/dis_run/saige/'
saige.files <- lapply(list.files(pattern = 'SAIGE_step2_chr.*.txt', 
                                 path = paste0(saige.dir, 'outputs'), 
                                 full.names = T), 
                      read.table, comment.char = "", header = T)
saige.names <- lapply(saige.files, function(x) paste0('chr', x$CHR[1], sep = ''))
names(saige.files) <- saige.names

saige <- Reduce(rbind, saige.files)
write.table(saige, paste0(saige.dir, 'saige_summary_statistics.txt'),
            row.names = F, quote = F, sep = '\t')
R.utils::gzip(paste0(saige.dir, 'saige_summary_statistics.txt'), destname = paste0(saige.dir, 'saige_summary_statistics.txt.gz'), 
              overwrite = F, remove = F)

saige.sig <- filter(saige, p.value < 5E-8)
write.table(saige.sig, paste0(saige.dir, 'saige_sig_summary_statistics.txt'),
            row.names = F, quote = F, sep = '\t')

saige.log10 <- saige %>%
  mutate(log10Praw = -1*log10(p.value)) %>%
  mutate(log10P = ifelse(log10Praw > 40, 40, log10Praw)) %>%
  mutate(log10Plevel = ifelse(p.value < 5E-08, "possible", NA)) %>%
  mutate(log10Plevel = ifelse(p.value < 5E-09, "likely", log10Plevel))

# Reduction of the GWAS object
# This is for more efficient plotting
# This drops everything not useful in future AUC calcs estiamte in Nalls et al., 2019

saige.filt <- filter(saige.log10, log10P > 3.114074) %>%
  arrange(CHR, POS)

saige.plotting <- saige.log10 %>%
  group_by(CHR) %>%
  summarise(chr.len=max(POS)) %>%
  mutate(tot = cumsum(as.numeric(chr.len)) - chr.len) %>%
  dplyr::select(-chr.len) %>%
  left_join(saige.filt, ., by = 'CHR') %>%
#  left_join(saige.log10, ., by = 'CHR') %>%
  arrange(ordered(CHR), POS) %>%
  mutate(POScum=POS+tot)

saige.axis.df <- saige.plotting %>%
  group_by(CHR) %>%
  summarise(center = (max(POScum) + min(POScum)) / 2)

saige.manhattan <-
  ggplot(saige.plotting, aes(POScum, log10P)) +
  geom_point(aes(colour=as.factor(CHR)), alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("dark grey", "black"), 22)) +
  scale_x_continuous(label = saige.axis.df$CHR, breaks = saige.axis.df$center) +
  scale_y_continuous(expand = c(0, 0.05)) +
  theme_classic() +
  labs(x = "CHR", y = "-log10P") +
  theme(legend.position = "none") +
  geom_abline(intercept = -log10(5E-08), slope = 0, linetype = 2)
ggsave(paste0(saige.dir, 'saige_manhattan.png'), plot = saige.manhattan, dpi = 400)

saige.qqplotter <- function(df) {
  df <- df[!is.na(df$p.value), ]
  df <- mutate(df, exp.pvalues = (rank(df$p.value, ties.method="first")+.5)/(length(df$p.value)+1))
  ggplot(df, aes(-log10(exp.pvalues), -log10(p.value))) + 
    geom_point() +
    labs(x = "Expected logP", y = "Observed logP") + 
    geom_abline(slope = 1, intercept = 0) +
    #geom_label_repel(data = arrange(df, p.value)[1:5, ], box.padding = 0.5) +
    theme_classic() 
}

saige.qqplt <- saige.qqplotter(saige)
# ggsave('C:/Users/Prasanth/Documents/cardiac_gwas/clean_gwas/saige_qqplot.png',
#        plot = plink.qqplt, dpi = 400)

saige.chisq <- qchisq(1-saige$p.value, 1)
saige.lambda <- median(saige.chisq) / qchisq(0.5,1)
write.table(saige.lambda, paste0(saige.dir, 'saige_lambda_value.txt'))

saige.alt.chisq <- qnorm(saige$p.value/2)
saige.alt.lambda <- median(saige.alt.chisq^2, na.rm=T)/qchisq(0.5,df=1)
