
rm(list=ls())

library(plyr)
library(ggplot2)
library(latex2exp)


# Load and organise data
loo.results = read.csv(file.path("data", "result_summaries", "SIM_loo_results.csv"))
loo.results$err = loo.results$err / 43
names(loo.results)[names(loo.results) == 'err'] = 'Eloo'

tru.results = read.csv(file.path("data", "result_summaries", "SIM_tru_results.csv"))
names(tru.results)[names(tru.results) == 'err'] = 'Etrue'
names(tru.results)[names(tru.results) == 'n.dims'] = 'n.dims.tru'

df = merge(loo.results, tru.results)
df[df$Etrue < 1e-15, 'Etrue'] = 1e-15

df.CV3 = subset(df, CV == 3)
df = subset(df, CV != 3)

df.CV2 = subset(df, CV == 2)

tmp1 = df.CV2
tmp1$n.dims = tmp1$n.dims.min
tmp1$range = 'min'
tmp1 = tmp1[, c('sim', 'method', 'Pi', 'Alt', 'n.dims', 'range')]

tmp2 = df.CV2
tmp2$n.dims = tmp2$n.dims.max
tmp2$range = 'max'
tmp2 = tmp2[, c('sim', 'method', 'Pi', 'Alt', 'n.dims', 'range')]

df.CV2 = rbind(tmp1, tmp2)


# Opti
df.opti = subset(df, Alt %in% 1:3)[, c('sim', 'method', 'Alt', 'n.dims.tru', 'Etrue')]
names(df.opti)[names(df.opti) == 'n.dims.tru'] = 'n.dims'
df.opti$Alt = as.character(df.opti$Alt)
df.opti[df.opti$Alt == "1", 'Alt'] = 'CV2'
df.opti[df.opti$Alt == "2", 'Alt'] = 'CV1Pi1'
df.opti[df.opti$Alt == "3", 'Alt'] = 'CV1Pi2'

# E_{true, k}
Etrue = read.csv(file.path("output", "result_summaries", "SIM_Etrue.csv"))
Etrue$err = Etrue$fn + Etrue$fp
Etrue[Etrue$err < 1e-15, 'err'] = 1e-15
Etrue = Etrue[, c('sim', 'method', 'n.dims', 'err')]

kstar = ddply(Etrue, 
              c('sim', 'method'), 
              summarise, 
              n.dims = min(n.dims[err == min(err)]))
tmp = transform(merge(Etrue, kstar), Alt = 'Opti')
names(tmp)[names(tmp) == 'err'] = 'Etrue'
df.opti = rbind(df.opti, tmp)



# Plotting aesthetic modifications
x_max_hist = max(c(max(df$Eloo*43), max(df$Etrue*43)))
df$CV = paste0('CV', df$CV)
df$Pi = paste0('Pi', df$Pi)
df.CV2$Pi = paste0('Pi', df.CV2$Pi)






# Figure 2: Ecv vs Epred
tmp = df[, c('sim', 'method', 'CV', 'Pi', 'Eloo')]
names(tmp)[names(tmp) == 'Eloo'] = 'E'
tmp$Type = 'CV'
df.plot = tmp

tmp = df[, c('sim', 'method', 'CV', 'Pi', 'Etrue')]
names(tmp)[names(tmp) == 'Etrue'] = 'E'
tmp$Type = 'Predict'
df.plot = rbind(df.plot, tmp)

tmp = df[, c('sim', 'method', 'CV', 'Pi', 'Eloo', 'Etrue')]
tmp$E = tmp$Eloo - tmp$Etrue
tmp$Type = 'CV - Predict'
df.plot = rbind(df.plot, tmp[, c('sim', 'method', 'CV', 'Pi', 'E', 'Type')])

df.plot$Alt = factor(paste(df.plot$CV, df.plot$Pi, sep = '.'), levels = c("CV1.Pi2", "CV1.Pi1", "CV2.Pi2", "CV2.Pi1"))

p = ggplot(df.plot, aes(x = E, fill = method)) +
  geom_histogram(binwidth = 1/43, 
                 position = 'identity', alpha = 0.4) +
  facet_grid(Alt ~ Type, scales = 'free')
ggsave(file.path('..', 'figures', 'Fig2_EcvEpred.png'), p, width = 7, height = 7)




# Figure 3 (Supplementary): kstar
tmp = df.CV2[, c('sim', 'method', 'Pi', 'range', 'n.dims')]
names(tmp)[names(tmp) == 'range'] = 'CV'
tmp$CV = paste0('CV2 (', tmp$CV, ')')
tmp = rbind(tmp, subset(df, CV == 'CV1')[, c('sim', 'method', 'Pi', 'CV', 'n.dims')])
df.plot = tmp

p = ggplot(df.plot, aes(x = n.dims, fill = method)) +
  geom_histogram(breaks = seq(1.5, 39.5, 2), position = 'identity', alpha = 0.4) +
  facet_grid(CV ~ Pi) + 
  xlab(TeX('$k^*$ used to calculate $E^{cv}$'))
# print(p)
ggsave(file.path('..', 'figures', 'Fig3_kstar.png'), p, width = 7, height = 7)






# Figure 4 (Supplementary): Etrue vs k
p = ggplot(Etrue, aes(x = n.dims, y = err, group = interaction(sim, method), colour = method)) +
  geom_line(alpha = 0.08, size = 0.6) + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) + 
  xlab('k') + 
  ylab(TeX('$E^{pred}_k$'))
# print(p)
ggsave(file.path('..', 'figures', 'Fig4_Epred_vs_k.png'), p, width = 7, height = 7)

