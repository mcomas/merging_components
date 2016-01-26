library(mixtools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(mixpack)

if(!file.exists('Pigs.rda')){
  download.file('https://github.com/cran/zCompositions/blob/master/data/Pigs.rda?raw=true',
                destfile='Pigs.rda')
}
load('Pigs.rda')

pigs <- Pigs %>% tbl_df

bic = function(loglik, k, n, d) -2 * loglik + (k * d - 1) * log(n)

set.seed(2)
K = 6
suppressMessages(
  FIT <- lapply(1:100, function(it, K){
    multmixEM(as.matrix(pigs), k = K)
  }, K)
)

fit = FIT[[which.min(sapply(FIT, function(fit) bic(fit$loglik, K, NROW(pigs), NCOL(pigs)) ))]]




d = (pigs.g <- pigs %>% mutate(
  TOTAL = BED+HALF.BED+PASSAGE+HALF.PASS+FEEDER+HALF.FEED,
  Comp = as.character(apply(fit$posterior, 1, which.max)))) %>%
  group_by(Comp) %>%
  summarise(
    BED = mean(BED/TOTAL),
    HALF.BED = mean(HALF.BED/TOTAL),
    PASSAGE = mean(PASSAGE/TOTAL),
    HALF.PASS = mean(HALF.PASS/TOTAL),
    FEEDER = mean(FEEDER/TOTAL),
    HALF.FEED = mean(HALF.FEED/TOTAL)) %>%
  gather(key=Variable, value=Proportion, -Comp)

ggplot() +
  geom_bar(data = d, aes(x=Variable, y=Proportion, fill=Comp), color='black', stat = 'identity') +
  facet_wrap(~Comp, ncol=3) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("")

ggsave(filename = 'figures/multinomial_mixt_all.pdf', height=6)




hp = get_hierarchical_partition(fit$posterior, omega = 'prop', lambda = 'coda.norm')
seq_partition = lapply(hp, sapply, paste, collapse=',')
cat(paste(sapply(seq_partition, function(partition) paste(sprintf('\\{%s\\}', partition), collapse=',')), collapse='\\}, \\\\ \n \\mathcal{P}_ &=& \\;\\; \\{'))

df = data.frame(
  Clusters = 1:6,
  S.values = attr(hp, 'S.value')
)
ggplot() +
  geom_point(data=df, aes(x=Clusters, y=S.values),size=2) + 
  theme_classic() +
  ggtitle('S-values during the merging process') +
  xlab('Clusters') + ylab('S-value')
ggsave(filename = 'figures/multinomial_Svalues_all.pdf', width = 7, height=4.5)
# pdf(file = 'figures/multinomial_Svalues.pdf', height=5.5)
# plot(attr(hp, 'S.value'), xlab='Clusters', ylab='S-value', main = 'S-values during the merging process')
# dev.off()

K = 3
clusters = mixpack::cluster_partition(fit$posterior, partition = hp[[K]] %>% setNames(., sapply(., function(lbl) paste(lbl, collapse='-'))))
d = pigs.g %>%
  mutate(
    Cluster = clusters
  ) %>% group_by(Cluster, Comp) %>%
  summarise(
    BED = mean(BED/TOTAL),
    HALF.BED = mean(HALF.BED/TOTAL),
    PASSAGE = mean(PASSAGE/TOTAL),
    HALF.PASS = mean(HALF.PASS/TOTAL),
    FEEDER = mean(FEEDER/TOTAL),
    HALF.FEED = mean(HALF.FEED/TOTAL))%>%
  gather(key=Variable, value=Proportion, -Cluster, -Comp)

ggplot() +
  geom_bar(data = d, aes(x=Variable, y=Proportion, fill=Comp), color='black', stat = 'identity', position=position_dodge()) +
  facet_wrap(~Cluster, ncol=3) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("")
ggsave(filename = 'figures/multinomial_clust3_all.pdf', height = 3.5)