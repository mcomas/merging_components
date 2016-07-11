f_omega1 = function(v_tau, a) as.numeric(which.max(v_tau) == a)
f_omega2 = function(v_tau, a) v_tau[a]
f_lambda= function(v_tau, a, b) exp(-log(v_tau[a]/v_tau[b])^2)


HP1a = get_hierarchical_partition(fit$posterior, omega = 'prop', lambda = 'demp')
HP2a = get_hierarchical_partition(fit$posterior, omega = 'prop', lambda = 'demp.mod')
HP3a = get_hierarchical_partition(fit$posterior, omega = 'cnst', lambda = 'entr')
HP4a = get_hierarchical_partition(fit$posterior, omega = 'prop', lambda = 'coda')


HP1b = get_hierarchical_partition(fit$posterior, omega = 'prop', lambda = 'prop')
HP2b = get_hierarchical_partition(fit$posterior, f_omega = f_omega2, f_lambda = f_lambda)





HP3a = get_hierarchical_partition(fit$posterior, omega = 'cnst', lambda = 'entr')
HP3b = get_hierarchical_partition(fit$posterior, f_omega = f_omega2, f_lambda = f_lambda)


HP5a = get_hierarchical_partition(fit$posterior, omega = 'prop', lambda = 'coda.norm')
HP6 = get_hierarchical_partition(fit$posterior, omega = 'prop', lambda = 'coda')

df2 = data_frame(
  K = factor(2:6, rev(2:6)),
  'DEMP' = attr(HP1a, 'S.value')[1:5],
  'DEMP modified' = attr(HP2a, 'S.value')[1:5],
  'Proportion' = attr(HP1b, 'S.value')[1:5],
  'Aitchison (scaled)' = 1-exp(-attr(HP2b, 'S.value')[1:5]),
  'Dif. Entropy (scaled)' = attr(HP3a,'S.value')[1:5]/(-log(1/2:6)),
  'Log-ratio (scaled)' = exp(attr(HP4a, 'S.value')[1:5])/(1 + exp(attr(HP4a, 'S.value')[1:5]))
) %>% gather(key=method, value=S.value, -K)
df2$method = factor(df2$method, levels = c('DEMP', 'DEMP modified', 'Proportion', 'Aitchison (scaled)', 
                                           'Dif. Entropy (scaled)', 'Log-ratio (scaled)'))

ggplot() +
  geom_point(data=df2, aes(x=K, y=S.value), size=3) +
  facet_wrap(~method, scales = 'free', nrow=2 ) +
  theme_classic() +
  xlab('Clusters') + ylab('S-value') +
  ggtitle('S-value during the merging process') +
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"))
ggsave(filename = 'figures/multinomial_Svalues_all.pdf', width=7.5, height=5)
