NDATA = 500
set.seed(1)

library(MixSim)
SEP = 40
ms = list(Pi = rep(1/6, 6),
     Mu = array(c(0,  0,SEP,SEP,SEP,   0,
                  0,  0,SEP,SEP,  0, SEP), c(6,2)),
     S  = array(c(50, 0,
                  0, 5,
                  5, 0,
                  0, 50), c(2,2,6)))
       
SIM = simdataset(n=NDATA, Pi=ms$Pi, Mu=ms$Mu, S=ms$S)

df = data.frame(SIM$X)

plot(df)
library(Rmixmod)
mixt = mixmodCluster(df, nbCluster=6, strategy=mixmodStrategy(nbTry = 1, seed = 1),
                     models = mixmodGaussianModel(listModels="Gaussian_pk_Lk_Ck"))
ord = c(1,6,4,2,5,3)
# P  = round(mixt@bestResult@parameters@proportions, 2)
P = ms$Pi
# M  = round(mixt@bestResult@parameters@mean, 2)
M = ms$Mu[ord,]
# S = lapply(mixt@bestResult@parameters@variance, round, 2)
S = lapply(ord, function(i) ms$S[,,i])

library(stringr)
library(Hmisc)
sink(file = 'tex/partition-example-pars.tex')
for(i in 1:4){
  cat('\\[\n')
  cat('\\begin{array}{l@{\\hskip 0.1in}l@{\\hskip 0.1in}c }\n')
  
  cat(sprintf("\\pi_%d = \\frac{1}{6}, & \\m\\mu_%d = \\left(%s\\right), & \\m\\Sigma_%d = \\left(\n%s\\right), \\\\ & &\\\\ \n",
              i, 
              i, paste(M[i,], collapse = ','),
              i, str_replace_all(latexTabular(S[[i]], headings = NULL, helvetica = FALSE), 'tabular', 'array')))
  
  cat('\\end{array}\n')
  cat('\\]\n')
}
cat('\\[\n')
cat('\\begin{array}{l@{\\hskip 0.1in}l@{\\hskip 0.1in}c }\n')
i = 5
cat(sprintf("\\pi_%d = \\frac{1}{6}, & \\m\\mu_%d = \\left(%s\\right), & \\m\\Sigma_%d = \\left(\n%s\\right) \\text{and} \\\\ & &\\\\ \n",
            i, 
            i, paste(M[i,], collapse = ','),
            i, str_replace_all(latexTabular(S[[i]], headings = NULL, helvetica = FALSE), 'tabular', 'array')))

cat('\\end{array}\n')
cat('\\]\n')

cat('\\[\n')
cat('\\begin{array}{l@{\\hskip 0.1in}l@{\\hskip 0.1in}c }\n')
i = 6
cat(sprintf("\\pi_%d = \\frac{1}{6}, & \\m\\mu_%d = \\left(%s\\right), & \\m\\Sigma_%d = \\left(\n%s\\right). \\\\ & &\\\\ \n",
            i, 
            i, paste(M[i,], collapse = ','),
            i, str_replace_all(latexTabular(S[[i]], headings = NULL, helvetica = FALSE), 'tabular', 'array')))

cat('\\end{array}\n')
cat('\\]\n')
sink()

library(ggplot2)
library(latex2exp)

xrange = round(1.1*range(df$X1))
yrange = round(1.1*range(df$X2))
(p <- ggplot() + 
  geom_point(data=df, aes(x=X1, y=X2), col='gray') + 
  geom_point(data=data.frame(M), aes(x=X1, y=X2), col = 'black', shape = 3, size=7) + 
  coord_fixed(xlim=xrange, ylim=yrange) + 
  theme_bw() +
  xlab(latex2exp('x_1')) + ylab(latex2exp('x_2')))



xlimits = seq(-30,130, 0.5)
ylimits = seq(-30,130, 0.5)
cm = expand.grid(X1 = xlimits, X2 = ylimits)
library(mixpack)
cm$z = dmixnorm_solution(cm, solution=mixt)

p.c6 <- p + stat_contour(data=cm, aes(x=X1, y=X2, z=z), col='blue')  +
  geom_point(data=data.frame(M), aes(x=X1, y=X2), col = 'black', shape = 3, size=7)
ggsave(p.c6, filename = 'figures/partition-example-mixture.pdf', width = 7, height=6)

library(dplyr)
library(mvtnorm)
CN = lapply(1:6, function(i){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = P[i] * dmvnorm(cn, mean = M[i,], sigma = S[[i]])
  cn$id = sprintf('{%s}',i)
  cn
}) %>% bind_rows

ii = apply( mixt@bestResult@proba, 1, which.max)
df$id = sapply(ii, function(i) sprintf('{%s}',paste(i, collapse=',')) )

p.all <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN, aes(x=X1, y=X2, z=z), alpha=0.6, col='blue') + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none") +
  scale_x_continuous(breaks=c(0, 20, 40)) + scale_y_continuous(breaks=c(0, 20, 40))
ggsave(p.all, filename = 'figures/partition-example-part6.pdf', width = 10, height=2.5)

library(tidyr)

HP = get_hierarchical_partition(mixt@bestResult@proba[,ord], omega = 'cnst', lambda = 'entr')
HP2 = get_hierarchical_partition(mixt@bestResult@proba[,ord], omega = 'dich', lambda = 'demp.mod')
df = data.frame(
  K = 1:6,
  'Entropy' = attr(HP,'S.value'),
  'DEMP' = attr(HP2, 'S.value')
) %>% gather(key=method, value=S.value, -K)

HP = lapply(HP, lapply, function(v) ord[v])
HP2 = lapply(HP2, lapply, function(v) ord[v])

ggplot() +
  geom_point(data=df, aes(x=K, y=S.value)) +
  facet_wrap(~method, scales = 'free') +
  theme_classic() +
  xlab('Clusters') + ylab('S-value') +
  ggtitle('S-value during the merging process')
ggsave(filename = 'figures/gaussian_Svalues.pdf', height = 3.25)

partition = HP[[3]]
names(partition) = sapply(partition, function(v) paste(sort(v), collapse=','))

CN3 = lapply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_solution(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
}) %>% bind_rows

prop_partition = function(prop, partition) lapply(partition, function(part, prop){
  if(length(part) == 1){
    return(prop[,part])
  }else{
    apply(prop[,part], 1, sum)
  }
}, prop) %>%data.frame %>% bind_cols

ii = apply( prop_partition(mixt@bestResult@proba %>% as.data.frame, partition), 1, which.max)
df$id = sapply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn3 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN3, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none") +
  scale_x_continuous(breaks=c(0, 20, 40)) + scale_y_continuous(breaks=c(0, 20, 40))
ggsave(p.cn3, filename = 'figures/partition-example-part3a.pdf', width = 5, height=2.5)

partition = list(c(1, 2, 3), 4, c(5, 6))

CN3b = lapply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_solution(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
}) %>% bind_rows

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = sapply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn3b <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN3b, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none") +
  scale_x_continuous(breaks=c(0, 20, 40)) + scale_y_continuous(breaks=c(0, 20, 40))
ggsave(p.cn3b, filename = 'figures/partition-example-part3b.pdf', width = 5, height=2.5)

partition = HP[[5]]

CN5 = lapply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_solution(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
}) %>% bind_rows

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = sapply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn5 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN5, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none") +
  scale_x_continuous(breaks=c(0, 20, 40)) + scale_y_continuous(breaks=c(0, 20, 40))
ggsave(p.cn5, filename = 'figures/partition-example-part5.pdf', width = 8.3, height=2.5)

partition = HP[[4]]

CN4 = lapply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_solution(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
}) %>% bind_rows

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = sapply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn4 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN4, aes(x=X1, y=X2, z=z), col='blue', alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none") +
  scale_x_continuous(breaks=c(0, 20, 40)) + scale_y_continuous(breaks=c(0, 20, 40))
ggsave(p.cn4, filename = 'figures/partition-example-part4.pdf', width = 6.6, height=2.5)

partition = HP[[2]]

CN2 = lapply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_solution(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
}) %>% bind_rows

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = sapply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn2 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN2, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none") +
  scale_x_continuous(breaks=c(0, 20, 40)) + scale_y_continuous(breaks=c(0, 20, 40))
ggsave(p.cn2, filename = 'figures/partition-example-part2.pdf', width = 3.3, height=2.5)

partition = HP[[1]]

CN1 = lapply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_solution(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
}) %>% bind_rows

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = sapply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn1 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN1, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none") +
  scale_x_continuous(breaks=c(0, 20, 40)) + scale_y_continuous(breaks=c(0, 20, 40))
ggsave(p.cn1, filename = 'figures/partition-example-part1.pdf', width = 2, height=2.5)
