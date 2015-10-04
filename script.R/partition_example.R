if(!exists('SEED')) SEED = 70
if(!exists('MAXOVERLAP')) MAXOVERLAP = 0.01
if(!exists('K0')) K0 = 3
if(!exists('DIM')) DIM = 2

if(!exists('NDATA')) NDATA = 500
if(!exists('Kf')) Kf = 6  ## Clusters utilitzats en el darrer pas

# if(!exists('NSIM')) NSIM = 100
# if(!exists('OFILE')) OFILE = 'data/sim-01.RData'

# library(plyr)
# library(reshape2)
# library(ggplot2)
# library(stringr)
# library(Hmisc)
# library(abind)

set.seed(SEED)

library(MixSim)
ms = MixSim(MaxOmega=MAXOVERLAP, K=K0, p=DIM, PiLow=0.01, sph=TRUE)
SIM = simdataset(n=NDATA, Pi=ms$Pi, Mu=ms$Mu, S=ms$S)

df = data.frame(100*SIM$X)

library(Rmixmod)
mixt = mixmodCluster(df, nbCluster=Kf, strategy=mixmodStrategy(nbTry = 10, seed = 1))

P  = round(mixt@bestResult@parameters@proportions, 2)
M  = round(mixt@bestResult@parameters@mean, 2)
S = lapply(mixt@bestResult@parameters@variance, round, 2)

library(stringr)
library(Hmisc)
sink(file = 'tex/partition-example-pars.tex')
for(i in 1:4){
  cat('\\[\n')
  cat('\\begin{array}{l@{\\hskip 0.1in}l@{\\hskip 0.1in}c }\n')
  
  cat(sprintf("\\hat{\\pi}_%d = %s, & \\hat{\\m\\mu}_%d = \\left(%s\\right), & \\hat{\\m\\Sigma}_%d = \\left(\n%s\\right), \\\\ & &\\\\ \n",
              i, paste(P[i], collapse = ','),
              i, paste(M[i,], collapse = ','),
              i, str_replace_all(latexTabular(S[[i]], headings = NULL, helvetica = FALSE), 'tabular', 'array')))
  
  cat('\\end{array}\n')
  cat('\\]\n')
}
cat('\\[\n')
cat('\\begin{array}{l@{\\hskip 0.1in}l@{\\hskip 0.1in}c }\n')
i = 5
cat(sprintf("\\hat{\\pi}_%d = %s, & \\hat{\\m\\mu}_%d = \\left(%s\\right), & \\hat{\\m\\Sigma}_%d = \\left(\n%s\\right) \\text{and} \\\\ & &\\\\ \n",
            i, paste(P[i], collapse = ','),
            i, paste(M[i,], collapse = ','),
            i, str_replace_all(latexTabular(S[[i]], headings = NULL, helvetica = FALSE), 'tabular', 'array')))

cat('\\end{array}\n')
cat('\\]\n')

cat('\\[\n')
cat('\\begin{array}{l@{\\hskip 0.1in}l@{\\hskip 0.1in}c }\n')
i = 6
cat(sprintf("\\hat{\\pi}_%d = %s, & \\hat{\\m\\mu}_%d = \\left(%s\\right), & \\hat{\\m\\Sigma}_%d = \\left(\n%s\\right). \\\\ & &\\\\ \n",
            i, paste(P[i], collapse = ','),
            i, paste(M[i,], collapse = ','),
            i, str_replace_all(latexTabular(S[[i]], headings = NULL, helvetica = FALSE), 'tabular', 'array')))

cat('\\end{array}\n')
cat('\\]\n')
sink()

library(ggplot2)
library(latex2exp)

#xrange = round(1.1*range(df$X1))
#yrange = round(1.1*range(df$X2))
ggplot() + 
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.3) + 
  geom_point(data=data.frame(M), aes(x=X1, y=X2), col = 'black', shape = 3, size=7) + 
  coord_fixed(xlim=xrange, ylim=yrange) + 
  theme_bw() +
  xlab(latex2exp('x_1')) + ylab(latex2exp('x_2'))



xlimits = seq(-30,130, 0.5)
ylimits = seq(-30,130, 0.5)
cm = expand.grid(X1 = xlimits, X2 = ylimits)
library(mixpack)
cm$z = dmixnorm_solution(cm, solution=mixt)

p.c6 <- p + stat_contour(data=cm, aes(x=X1, y=X2, z=z), col='blue')  +
  geom_point(data=data.frame(M), aes(x=X1, y=X2), col = 'black', shape = 3, size=7)
ggsave(p.c6, filename = 'figures/partition-example-mixture.pdf', width = 7, height=6)

CN = ldply(1:6, function(i){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = P[i] * dmvnorm(cn, mean = M[i,], sigma = S[[i]])
  cn$id = sprintf('{%s}',i)
  cn
})

ii = apply( mixt@bestResult@proba, 1, which.max)
df$id = laply(ii, function(i) sprintf('{%s}',paste(i, collapse=',')) )

p.all <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN, aes(x=X1, y=X2, z=z), alpha=0.6) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none") 
ggsave(p.all, filename = 'figures/partition-example-part6.pdf', width = 10, height=2.5)

partition = list(c(1, 3, 6), c(2, 5), 4)

CN3 = ldply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_rmixmod(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
})

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = laply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn3 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN3, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none")
ggsave(p.cn3, filename = 'figures/partition-example-part3a.pdf', width = 5, height=2.5)

partition = list(c(1, 2, 3), 4, c(5, 6))

CN3b = ldply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_rmixmod(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
})

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = laply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn3b <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN3b, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none")
ggsave(p.cn3b, filename = 'figures/partition-example-part3b.pdf', width = 5, height=2.5)

partition = list(c(1, 6), 2, 3, 4, 5)

CN5 = ldply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_rmixmod(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
})

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = laply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn5 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN5, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none")
ggsave(p.cn5, filename = 'figures/partition-example-part5.pdf', width = 8.3, height=2.5)

partition = list(c(1, 3, 6), 2, 4, 5)

CN4 = ldply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_rmixmod(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
})

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = laply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn4 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN4, aes(x=X1, y=X2, col=id, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none")
ggsave(p.cn4, filename = 'figures/partition-example-part4.pdf', width = 6.6, height=2.5)

partition = list(c(1, 3, 6), c(2, 4, 5))

CN2 = ldply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_rmixmod(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
})

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = laply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn2 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN2, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none")
ggsave(p.cn2, filename = 'figures/partition-example-part2.pdf', width = 3.3, height=2.5)

partition = list(c(1, 2, 3, 4, 5, 6))

CN1 = ldply(partition, function(part){
  cn = expand.grid(X1 = xlimits, X2 = ylimits)
  cn$z = dmixnorm_rmixmod(cn, mixt, part = part) #dmvnorm(cn, mean = Mu[,i], sigma = S[,,i])
  cn$id = sprintf('{%s}',paste(part, collapse=','))
  cn
})

ii = apply( prop_partition(mixt@bestResult@proba, partition), 1, which.max)
df$id = laply(ii, function(i) sprintf('{%s}',paste(partition[[i]], collapse=',')) )

p.cn1 <- ggplot() + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), alpha=0.1) +
  geom_point(data=df, aes(x=X1, y=X2), alpha=0.8, size=1) +
  stat_contour(data=CN1, aes(x=X1, y=X2, z=z), alpha=0.8) + 
  facet_wrap(~id, nrow=1) + theme_bw() + theme(legend.position="none")
ggsave(p.cn1, filename = 'figures/partition-example-part1.pdf', width = 2, height=2.5)
