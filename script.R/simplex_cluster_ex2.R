NDATA = 400
set.seed(1)
library(ggplot2)
library(MixSim)
SEP = 30
ms = list(Pi = c(2,1,1,2)/6,
          Mu = array(c(0+SEP/15,  SEP/2,SEP/2,SEP-SEP/15,
                       0,  SEP*sqrt(3)/2,SEP*sqrt(3)/2,  0), c(4,2)),
          S  = array(c(20, 0,
                       0, 20,
                       ###
                       40, 0,
                       0, 5,
                       5, 0,
                       0, 75,
                       ###
                       20, 0,
                       0, 20), c(2,2,4)))

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
cm$z = dmixnorm_solution(cm, solution=mixt)


library(tidyr)
step_merging = function(post, omega, lambda, f_omega = NULL, f_lambda = NULL){
  S.values = merge_step(post, omega, lambda, f_omega, f_lambda)
  
  S.values = S.values + diag(-Inf, NROW(S.values))
  ind = which(S.values == max(S.values), arr.ind = TRUE)
  new_post = merge_components(post, ind[1], ind[2])
  list('S.values' = S.values, 'post' = new_post)
}
L = list()
L[[1]] = step_merging(mixt@bestResult@proba[,ord], omega = 'prop', lambda = 'demp')
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
    Cluster = sprintf('g%d',apply(L[[1]]$post, 1, which.max)))
POST4 = POST4 %>%
  mutate(
    Cluster = sprintf('g%d',apply(mixt@bestResult@proba[,ord], 1, which.max)))
POST = POST %>%
  mutate(
    Cluster = sprintf('g%d',apply(L[[1]]$post, 1, which.max)))
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
p1 <- ggplot() + 
  geom_point(data=df, aes(x=X1, y=X2, col=Cluster), size=1.5) + 
  stat_contour(data=cm, aes(x=X1, y=X2, z=z), col='black', size=1) +
  coord_fixed(xlim=xrange, ylim=yrange) + 
  theme_bw() +
  xlab(TeX('x_1')) + ylab(TeX('x_2')) + 
  ggtitle('Sample') + ggplot.theme

p2a <- ggtern() +
  geom_point(data=POST4, aes(x = V1, y = V2, z = V3, col=Cluster), size=2.5) +
  ggtitle('Posterior probabilities') + 
  Tlab(expression('tau[2]')) +
  Llab(expression('tau[1]')) +
  Rlab(expression('tau[3]')) + ggtern.theme
p2b <- ggtern() +
  geom_point(data=POST4, aes(x = V1, y = V2, z = V4, col=Cluster), size=2.5) +
  ggtitle('Posterior probabilities') + 
  Tlab(expression('tau[2]')) +
  Llab(expression('tau[1]')) +
  Rlab(expression('tau[4]')) + ggtern.theme

options("tern.discard.external" = FALSE) 
p2 <- ggtern() +
  geom_point(data=POST, aes(x = V1, y = V2, z = V3, col=Cluster), size=3) +
  ggtitle('Posterior probabilities') + 
  Tlab(expression('tau[2]')) +
  Llab(expression('tau[1]')) +
  Rlab(expression('tau[3]')) + ggtern.theme
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

write.table(POST4, file = 'post4.txt', row.names=FALSE, sep=';')
write.table(POST, file='post3.txt', row.names=FALSE, sep=';')


p2bw <- ggtern() +
  geom_mask() +
  geom_point(data=POST, aes(x = V1, y = V2, z = V3), size=3) +
  ggtitle('Posterior probabilities') + 
  Tlab(expression('tau[2]')) +
  Llab(expression('tau[1]')) +
  Rlab(expression('tau[3]')) + ggtern.theme
p2bw
