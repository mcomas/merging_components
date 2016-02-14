library(fpc)
library(mixpack)
library(compositions)
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

c.stats.par = function(post, par){
  if(NCOL(post) == 1){
    return(NA)
  }
  post[post == 0] = .Machine$double.xmin
  apost = acomp(post)
  Adist = dist(apost)
  if(par == 'g2'){
    c.stats = cluster.stats(Adist, apply(post, 1, which.max), G2 = T, G3 = F)
  }else{
    c.stats = cluster.stats(Adist, apply(post, 1, which.max), G2 = F, G3 = F)
  }
  
  c.stats[[par]]
}

L = L.post(post = post, omega = 'prop', lambda = 'coda.norm')

L.index = function(L, index){
  sapply(L, function(l) c.stats.par(l$post, index))
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
V.index = lapply(names(INDEXS), function(index) L.index(L, index)) %>% as.data.frame %>% tbl_df
names(V.index) = INDEXS

df = V.index  %>% na.omit %>% mutate(
  K = factor(K, NROW(V.index):1)
) %>% gather(key = index, value = value, -K)

ggplot() +
  geom_line(data=df, aes(x=K, y=value), group=1) +
  geom_point(data=df, aes(x=K, y=value), size=3) +
  facet_wrap(~index, scales='free') + 
  theme_classic() +
  ylab('Statistic value') + xlab('Clusters')
ggsave(filename = 'figures/multinomial_statistics.pdf', width=9.5, height=7)


