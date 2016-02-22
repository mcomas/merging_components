library(dplyr)

library(ggplot2)
xlog = function(x) ifelse(x == 0, 0, x * log(x))

epsilon = 0.01
df.entropy = lapply(seq(0.2, 0.8, length.out=5), function(other){
  data_frame(
    tau_a = seq(epsilon, 1-other-epsilon, 0.0001),
    tau_b = 1-other-tau_a,
    index = xlog(tau_a+tau_b) - xlog(tau_a) - xlog(tau_b),
    other = as.character(other),
    lambda = 'Entropy'
  )
}) %>% bind_rows
df.aitchison = lapply(seq(0.2, 0.8, length.out=5), function(other){
  data_frame(
    tau_a = seq(epsilon/other, 1-other-epsilon/other, 0.0001),
    tau_b = 1-other-tau_a,
    index = -log(tau_b/tau_a)^2, #exp(-log(tau_b/tau_a)^2), #
    other = as.character(other),
    lambda = 'Aitchison distance'
  )
}) %>% bind_rows %>% subset(index > -5)
df = bind_rows(df.aitchison, df.entropy)

ggplot() +
  geom_line(data=df, aes(x=tau_b, y=index, linetype=other), size=1) +
  facet_wrap(~lambda, nrow=1, scale='free') +
  theme_classic() +
  theme(legend.position='top') +
  xlab('tau_b (or 1-tau_a)')
ggsave(file='figures/entr_dist.pdf', width = 9.5, height = 4.4)

df.demp2 = lapply(seq(0.2, 0.8, length.out=4), function(other){
  data_frame(
    tau_a = seq(epsilon/other, 1-other-epsilon/other, 0.0001),
    tau_b = 1-other-tau_a,
    index = tau_b/(tau_a+tau_b),
    other = as.character(other),
    lambda = 'DEMP 2'
  )
}) %>% bind_rows 
df.log = lapply(seq(0.2, 0.8, length.out=4), function(other){
  data_frame(
    tau_a = seq(epsilon/other, 1-other-epsilon/other, 0.0001),
    tau_b = 1-other-tau_a,
    index = log(tau_b/tau_a),
    other = as.character(other),
    lambda = 'log'
  )
})%>% bind_rows
df = bind_rows(df.demp2, df.log)

ggplot() +
  geom_line(data=df, aes(x=tau_b, y=index, linetype=other), size=1) +
  facet_wrap(~lambda, nrow=1, scale='free') +
  theme_classic() +
  theme(legend.position='top') +
  xlab('tau_b (or 1-tau_a)')
ggsave(file='figures/demp2_log.pdf', width = 9.5, height = 4.4)

df.demp = lapply(c(0.5, 0.3, 0.1), function(other){
  data_frame(
    tau_a = seq(0, 1-other, 0.0001),
    tau_b = 1-other-tau_a,
    index = as.numeric(tau_b > tau_a & tau_b > other),
    other = as.character(other),
    lambda = 'DEMP'
  )
}) %>% bind_rows
df.demp2 = lapply(c(0.5, 0.3, 0.1), function(other){
  data_frame(
    tau_a = seq(0, 1-other, 0.0001),
    tau_b = 1-other-tau_a,
    index = tau_b/(tau_a+tau_b),
    other = as.character(other),
    lambda = 'DEMP 2'
  )
}) %>% bind_rows 

df = bind_rows(df.demp, df.demp2)
ggplot() +
  geom_point(data=df, aes(x=tau_b, y=index, col=lambda), size=1) +
  facet_wrap(~other, nrow=1) +
  theme_classic() +
  xlab('tau_b (or 1-tau_a)')
