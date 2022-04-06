library(tidyverse)
library(ggrepel)

# read in fuma and magma files 
fuma.dir <- 'C:/Users/Prasanth/Documents/cardiac_gwas/dis_run/saige/fuma/'
fuma.files <- list.files(fuma.dir, pattern = '*.txt')
magma.files <- list.files(fuma.dir, pattern = '*.out')

readFuma <- function(file){
  readr::read_delim(paste0(fuma.dir, file), comment = '#', delim = '\t')
}

fuma <- lapply(fuma.files, readFuma)
names(fuma) <- fuma.files
magma <- lapply(magma.files, readFuma)
names(magma) <- magma.files

#############################################################################
# magma gene analysis

magma.plotting <- magma$magma.genes.out %>%
  group_by(CHR) %>%
  summarise(chr.len=max(START)) %>%
  mutate(tot = cumsum(as.numeric(chr.len)) - chr.len) %>%
  dplyr::select(-chr.len) %>%
  left_join(magma$magma.genes.out, ., by = 'CHR') %>%
  arrange(ordered(CHR), START) %>%
  mutate(POScum = START+tot, 
         sig = ifelse(P < 0.05/nrow(magma$magma.genes.out), 'SIG', 'NOT_SIG'))

magma.axis.df <- magma.plotting %>%
  group_by(CHR) %>%
  summarise(center = (max(POScum) + min(POScum)) / 2)

magma.manhattan <-
  ggplot(magma.plotting, aes(POScum, -log10(P))) +
  geom_point(aes(colour=as.factor(CHR)), alpha = 0.8, size = 0.8) +
  scale_color_manual(values = rep(c("dark grey", "black"), 22)) +
  scale_x_continuous(label = magma.axis.df$CHR, breaks = magma.axis.df$center) +
  scale_y_continuous(expand = c(0, 0.05)) +
  theme_classic() +
  labs(x = "CHR", y = "-log10P") +
  theme(legend.position = "none") +
  geom_abline(intercept = -log10(0.05/nrow(magma.plotting)), slope = 0, linetype = 2) +
  geom_label_repel(data = filter(magma.plotting, sig == 'SIG'), 
                  aes(label = SYMBOL), max.overlaps = Inf, size = 2, label.size = 0.1, 
                  min.segment.length = 0, seed = 1)
#ggsave(paste0(fuma.dir, '/magma_manhattan.png'), plot = magma.manhattan, dpi = 400)

magma.qqplotter <- function(df) {
  df <- df[!is.na(df$P), ]
  df <- mutate(df, exp.pvalues = (rank(df$P, ties.method="first")+.5)/(length(df$P)+1))
  ggplot(df, aes(-log10(exp.pvalues), -log10(P))) + 
    geom_point() +
    labs(x = "Expected logP", y = "Observed logP") + 
    geom_abline(slope = 1, intercept = 0) +
    #geom_label_repel(data = arrange(df, p.value)[1:5, ], box.padding = 0.5) +
    theme_classic() 
}

magma.qqplt <- magma.qqplotter(magma$magma.genes.out)
#ggsave(paste0(fuma.dir, '/magma_qqplot.png'), plot = magma.qqplt, dpi = 400)

magma.chisq <- qchisq(1-magma$magma.genes.out$P, 1)
magma.lambda <- median(magma.chisq) / qchisq(0.5,1)
#write.table(magma.lambda, paste0(fuma.dir, '/magma_lambda_value.txt'))


#####################################################################################
# gwas catalogue comparisons

View(fuma$gwascatalog.txt)
