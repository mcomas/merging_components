library(fpc)
library(mixpack)
library(compositions)
# post = fit$posterior
# post = mixt@bestResult@proba
## Merging
step_merging = function(post, omega, lambda, f_omega = NULL, f_lambda = NULL){
  S.values = merge_step(post, omega, lambda, f_omega, f_lambda)
  
  S.values = S.values + diag(-Inf, NROW(S.values))
  ind = which(S.values == max(S.values), arr.ind = TRUE)
  new_post = merge_components(post, ind[1], ind[2])
  list('S.values' = S.values, 'post' = new_post)
}

L.post = function(post, omega, lambda){
  L = list()
  L[[1]] = list('S.value' = NA, 'post' = post)
  for(i in 2:NCOL(post)){
    L[[i]] = step_merging(L[[i-1]]$post, omega, lambda)
  }
  L
}

c.stats.par = function(post, par, post0 = NULL){
  if(is.null(post0)){
    post0 = post
  }
  if(NCOL(post) == 1){
    return(NA)
  }
  post0[post0 == 0] = .Machine$double.xmin
  apost = acomp(post0)
  Adist = dist(apost)
  if(par == 'g2'){
    c.stats = cluster.stats(Adist, apply(post, 1, which.max), G2 = T, G3 = F)
  }else{
    c.stats = cluster.stats(Adist, apply(post, 1, which.max), G2 = F, G3 = F)
  }
  
  c.stats[[par]]
}

L = L.post(post = post, omega = 'prop', lambda = 'coda.norm')

L.index = function(L, index, post0){
  sapply(L, function(l) c.stats.par(l$post, index, post0 = post0))
}

INDEXS = c('cluster.number' = 'K', 
           'average.between' = 'Average Between', 
           'average.within' = 'Average Within', 
           'ch' = 'Calinski-Harabasz', 
#            'pearsongamma' = 'Normalized gamma', 
#            'dunn' = 'Dunn', 
#            'dunn2' = 'Dunn2', 
#            'avg.silwidth' = 'Average silouhette', 
#            'sindex' = 'Separation',
           'g2' = 'G2')
V.index2 = lapply(names(INDEXS), function(index) L.index(L, index, L[[1]]$post)) %>% as.data.frame %>% tbl_df
names(V.index2) = INDEXS

V.index = lapply(names(INDEXS), function(index) L.index(L, index, NULL)) %>% as.data.frame %>% tbl_df
names(V.index) = INDEXS

df2 = V.index2  %>% na.omit %>% mutate(
  K = factor(K, NROW(V.index2):1)
) %>% gather(key = index, value = value, -K)

df = V.index  %>% na.omit %>% mutate(
  K = factor(K, NROW(V.index):1)
) %>% gather(key = index, value = value, -K)

ggplot() +
  #geom_line(data=df, aes(x=K, y=value), group=1) +
  geom_point(data=df, aes(x=K, y=value), size=3) +
  facet_wrap(~index, scales='free', nrow=1) + 
  theme_classic() +
  ylab('Statistic value') + xlab('Clusters') +
  ggtitle('Posteriori probability based statistics')
#ggsave(filename = 'figures/multinomial_statistics.pdf', width=10, height=2.7)
#ggsave(filename = 'figures/gaussian_statistics.pdf', width=10, height=2.7)
ggplot() +
  #geom_line(data=df, aes(x=K, y=value), group=1) +
  geom_point(data=df2, aes(x=K, y=value), size=3) +
  facet_wrap(~index, scales='free', nrow=1) + 
  theme_classic() +
  ylab('Statistic value') + xlab('Clusters') +
  ggtitle('Posteriori probability based statistics (k=6)')
#ggsave(filename = 'figures/multinomial_statistics2.pdf', width=10, height=2.7)
