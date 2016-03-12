NDATA = 500
set.seed(1)
library(ggplot2)
library(MixSim)
SEP = 20
ms = list(Pi = c(1,1,1,1.5)/4.5,
          Mu = array(c(0,SEP,SEP/2,SEP/2,
                       0,  0,SEP*sqrt(3)/2,SEP*sqrt(3)/2), c(4,2)),
          S  = array(c(10, 0,
                       0, 10,
                       10, 0,
                       0, 10,
                       10, 0,
                       0, 10,
                       1, 0,
                       0, 50), c(2,2,4)))

SIM = simdataset(n=NDATA, Pi=ms$Pi, Mu=ms$Mu, S=ms$S)

df = data.frame(SIM$X)

plot(df)

library(Rmixmod)
mixt = mixmodCluster(df, nbCluster=4, strategy=mixmodStrategy(nbTry = 1, seed = 1),
                     models = mixmodGaussianModel(listModels="Gaussian_pk_Lk_Ck"))

ord = 1:4
# P  = round(mixt@bestResult@parameters@proportions, 2)
P = ms$Pi
# M  = round(mixt@bestResult@parameters@mean, 2)
M = ms$Mu[ord,]
# S = lapply(mixt@bestResult@parameters@variance, round, 2)
S = lapply(ord, function(i) ms$S[,,i])

library(ggplot2)
library(latex2exp)

xrange = round(1.1*range(df$X1))
yrange = round(1.1*range(df$X2))


xlimits = seq(-30,130, 0.5)
ylimits = seq(-30,130, 0.5)
cm = expand.grid(X1 = xlimits, X2 = ylimits)
library(mixpack)
cm$z = dmixnorm(cm[,1:2], Pi = Pi, Mu = t(ms$Mu), S = ms$S)
cm$z1 = Pi[1] * mvtnorm::dmvnorm(cm[,1:2], mean = ms$Mu[1,], sigma = ms$S[,,1])
cm$z2 = Pi[2] * mvtnorm::dmvnorm(cm[,1:2], mean = ms$Mu[2,], sigma = ms$S[,,2])
cm$z3 = Pi[3] * mvtnorm::dmvnorm(cm[,1:2], mean = ms$Mu[3,], sigma = ms$S[,,3])
cm$z4 = Pi[4] * mvtnorm::dmvnorm(cm[,1:2], mean = ms$Mu[4,], sigma = ms$S[,,4])
cm$z34 = sum(Pi[3:4]) * dmixnorm(cm[,1:2], Pi = Pi[3:4], Mu = t(ms$Mu)[,3:4], S = ms$S[,,3:4])

library(tidyr)
step_merging = function(post, omega, lambda, f_omega = NULL, f_lambda = NULL){
  S.values = merge_step(post, omega, lambda, f_omega, f_lambda)
  
  S.values = S.values + diag(-Inf, NROW(S.values))
  ind = which(S.values == max(S.values), arr.ind = TRUE)
  new_post = merge_components(post, ind[1], ind[2])
  list('S.values' = S.values, 'post' = new_post)
}
L = list()
L[[1]] = step_merging(mixt@bestResult@proba[,ord], omega = 'prop', lambda = 'coda.norm')

library(dplyr)
POST4 = as.data.frame(mixt@bestResult@proba) %>% tbl_df
POST = as.data.frame(L[[1]]$post) %>% tbl_df
library(compositions)
aPOST = as.data.frame(scale(acomp(POST)))
library(ggtern)
ilrPOST = mixpack::ilr_coordinates(POST)
library(gtable)
library(latex2exp)
df = df %>%
  mutate(
    Cluster = sprintf('g%d',apply(POST, 1, which.max)),
    Cluster4 = sprintf('g%d',apply(POST4, 1, which.max)))
POST = POST %>%
  mutate(
    Cluster = sprintf('g%d',apply(L[[1]]$post, 1, which.max)))
POST4 = POST4 %>%
  mutate(
    Cluster = sprintf('g%d',apply(mixt@bestResult@proba, 1, which.max)))
aPOST = aPOST %>%
  mutate(
    Cluster = sprintf('g%d',apply(L[[1]]$post, 1, which.max)))
ilrPOST = ilrPOST %>%
  mutate(
    Cluster = sprintf('g%d',apply(L[[1]]$post, 1, which.max)))

ggplot.theme = ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title = element_text(size=17),
    legend.position = 'none')
ggtern.theme = ggtern::theme_bw() + 
  ggtern::theme(
    tern.axis.title = element_text(size=17),
    legend.position = 'none')
p1a <- ggplot() + 
  geom_point(data=df, aes(x=X1, y=X2, fill=Cluster4, shape=Cluster4), colour='black', size=2.2) +
  #geom_point(data=df, aes(x=X1, y=X2, shape=Cluster, col=Cluster), size=1.7) +
  stat_contour(data=cm, aes(x=X1, y=X2, z=z4, col='g2'), linejoin='mitre', alpha=0.9, size=0.5, linetype=5) +
  stat_contour(data=cm, aes(x=X1, y=X2, z=z3, col='g4'), linejoin='mitre', alpha=0.9, size=0.5, linetype=5) +
  stat_contour(data=cm, aes(x=X1, y=X2, z=z1, col='g1'), linejoin='mitre', alpha=0.9, size=0.5, linetype=5) +
  stat_contour(data=cm, aes(x=X1, y=X2, z=z2, col='g3'), linejoin='mitre', alpha=0.9, size=0.5, linetype=5) +
  #stat_contour(data=cm, aes(x=X1, y=X2, z=z), col='black', linejoin='mitre', alpha=0.8, size=0.5, linetype=5) +
  coord_fixed(xlim=xrange, ylim=yrange) + 
  scale_shape_manual(values = c(21,22,24,25)) +
  theme_bw() +
  xlab(TeX('x_1')) + ylab(TeX('x_2')) + 
  ggtitle('Sample') + ggplot.theme
ggsave(p1a, filename = 'figures/ex_4clust.pdf', width = 4, height=4)

p1 <- ggplot() + 
  geom_point(data=df, aes(x=X1, y=X2, fill=Cluster, shape=Cluster), colour='black', size=2.2) +
  #geom_point(data=df, aes(x=X1, y=X2, shape=Cluster, col=Cluster), size=1.7) +
  stat_contour(data=cm, aes(x=X1, y=X2, z=z34, col='g2'), linejoin='mitre', alpha=0.9, size=0.5, linetype=5) +
  stat_contour(data=cm, aes(x=X1, y=X2, z=z1, col='g1'), linejoin='mitre', alpha=0.9, size=0.5, linetype=5) +
  stat_contour(data=cm, aes(x=X1, y=X2, z=z2, col='g3'), linejoin='mitre', alpha=0.9, size=0.5, linetype=5) +
  #stat_contour(data=cm, aes(x=X1, y=X2, z=z), col='black', linejoin='mitre', alpha=0.8, size=0.5, linetype=5) +
  coord_fixed(xlim=xrange, ylim=yrange) + 
  scale_shape_manual(values = c(21,22,24)) +
  theme_bw() +
  xlab(TeX('x_1')) + ylab(TeX('x_2')) + 
  ggtitle('Sample') + ggplot.theme
ggsave(p1, filename = 'figures/ex_3clust.pdf', width = 4, height=4)
p2 <- ggtern() +
  geom_mask() +
  geom_point(data=POST, aes(x = V1, y = V2, z = V3, fill=Cluster, shape=Cluster), colour='black', size=2.2) +
  scale_shape_manual(values = c(21,22,24)) +
  #geom_point(data=POST, aes(x = V1, y = V2, z = V3, shape=Cluster, col=Cluster), size=1.7) +
  ggtitle('Posterior probabilities') + 
  Tlab(expression('tau[2]')) +
  Llab(expression('tau[1]')) +
  Rlab(expression('tau[3]')) + ggtern.theme
ggsave(p2, filename = 'figures/ex_ternary.pdf', width = 3.5, height=3.5)

p2b <- ggtern() +
  geom_mask() +
  geom_point(data=POST, aes(x = V1, y = V2, z = V3, fill=Cluster, shape=Cluster), colour='black', size=2.2) +
  scale_shape_manual(values = c(21,22,24)) +
  #geom_point(data=POST, aes(x = V1, y = V2, z = V3, shape=Cluster, col=Cluster), size=1.7) +
  ggtitle('Posterior probabilities') + 
  Tlab(expression('tau[2]')) +
  Llab(expression('tau[1]')) +
  Rlab(expression('tau[4]')) + ggtern.theme
ggsave(p2b, filename = 'tempo.pdf', width = 3.5, height=3.5)
p3 <- ggtern() +
  geom_point(data=aPOST, aes(x = V1, y = V2, z = V3, col=Cluster)) +
  ggtitle('Posterior probabilities after scaling') + 
  Tlab(expression('tau[2]^c')) +
  Llab(expression('tau[1]^c')) +
  Rlab(expression('tau[3]^c')) + ggtern.theme
p4 <- ggplot() +
  geom_point(data = ilrPOST, aes(x = coord.1, y = coord.2, col=Cluster)) +
  ggtitle('ILR coordinates') + ggplot2::theme_bw() +
  xlab( expression(paste( sqrt(1 / 2),'  ', ln(paste(' ',  tau[1] / tau[3], ' ') )) )) + 
  ylab( expression(paste( sqrt(2 / 3),'  ', ln(paste(' ',  sqrt(paste(tau[1],'·',  tau[3])) / tau[2],' ') )) )) +
  ggplot.theme



library(grid)
grid.newpage() 
vpa_ <- viewport(width = 0.5, height = 0.5, x = 0.25, y = 0.75) 
vpb_ <- viewport(width = 0.5, height = 0.5, x = 0.75, y = 0.75) 
vpc_ <- viewport(width = 0.5, height = 0.5, x = 0.25, y = 0.25) 
vpd_ <- viewport(width = 0.5, height = 0.5, x = 0.75, y = 0.25) 
print(p1, vp = vpa_) 
print(p2, vp = vpb_) 
print(p3, vp = vpc_) 
print(p4, vp = vpd_)

POST3  = POST %>% arrange(Cluster) %>% data.frame
POST4  = POST4 %>% arrange(Cluster) %>% data.frame
write.table(POST4, file='posterior4.csv', row.names = FALSE, sep = ';')
WriteXLS::WriteXLS(c('POST3', 'POST4'), ExcelFileName = 'posteriors.xls')
