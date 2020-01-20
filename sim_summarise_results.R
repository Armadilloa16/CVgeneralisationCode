
library(plyr)

N = 200

# For each simlation, load the results.
clas.pop = data.frame()
clas.tt = data.frame()
loo.alt1.cv = data.frame()
loo.alt1.tt = data.frame()
loo.alt2.cv = data.frame()
loo.alt2.tt = data.frame()
loo.cv6 = data.frame()
for (sim in 1:N) {
  cat(sim, ' ')
  
  clas.pop = rbind(clas.pop,
                   transform(read.csv(file.path("data", "sim_results",
                                                paste0("cca_lda_pop_",
                                                       toString(sim),
                                                       ".csv"))),
                             sim = sim, 
                             method = 'CCA-LDA'))
  clas.pop = rbind(clas.pop,
                   transform(read.csv(file.path("data", "sim_results",
                                                paste0("pca_lda_pop_",
                                                       toString(sim),
                                                       ".csv"))),
                             sim = sim,
                             method = 'PCA-LDA'))



  clas.tt = rbind(clas.tt,
                  transform(read.csv(file.path("data", "sim_results",
                                               paste0("cca_lda_tt_",
                                                      toString(sim),
                                                      ".csv"))),
                            sim = sim, 
                            method = 'CCA-LDA'))
  clas.tt = rbind(clas.tt,
                  transform(read.csv(file.path("data", "sim_results",
                                               paste0("pca_lda_tt_",
                                                      toString(sim),
                                                      ".csv"))),
                            sim = sim,
                            method = 'PCA-LDA'))



  loo.alt1.cv = rbind(loo.alt1.cv,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("cca_lda_alt1_loo_cv_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim, 
                                method = 'CCA-LDA'))
  loo.alt1.cv = rbind(loo.alt1.cv,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("pca_lda_alt1_loo_cv_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim,
                                method = 'PCA-LDA'))



  loo.alt1.tt = rbind(loo.alt1.tt,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("cca_lda_alt1_loo_tt_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim, 
                                method = 'CCA-LDA'))
  loo.alt1.tt = rbind(loo.alt1.tt,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("pca_lda_alt1_loo_tt_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim,
                                method = 'PCA-LDA'))



  loo.alt2.cv = rbind(loo.alt2.cv,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("cca_lda_alt2_loo_cv_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim, 
                                method = 'CCA-LDA'))
  loo.alt2.cv = rbind(loo.alt2.cv,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("pca_lda_alt2_loo_cv_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim,
                                method = 'PCA-LDA'))



  loo.alt2.tt = rbind(loo.alt2.tt,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("cca_lda_alt2_loo_tt_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim, 
                                method = 'CCA-LDA'))
  loo.alt2.tt = rbind(loo.alt2.tt,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("pca_lda_alt2_loo_tt_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim,
                                method = 'PCA-LDA'))


  
  
  
  loo.cv6 = rbind(loo.cv6,
                  transform(read.csv(file.path("data", "sim_results",
                                               paste0("cca_lda_cv6_loo_cv_",
                                                      toString(sim),
                                                      ".csv"))),
                            sim = sim, 
                            method = 'CCA-LDA'))
  loo.cv6 = rbind(loo.cv6,
                  transform(read.csv(file.path("data", "sim_results",
                                               paste0("pca_lda_cv6_loo_cv_",
                                                      toString(sim),
                                                      ".csv"))),
                            sim = sim,
                            method = 'PCA-LDA'))
  
  
}

# E_{true, k}
prior = 16/43
clas.pop = transform(clas.pop,
                  fn = prior       * pnorm(cutoff, mean = mu.p, sd = sd.p, lower.tail = TRUE),
                  fp = (1 - prior) * pnorm(cutoff, mean = mu.n, sd = sd.n, lower.tail = FALSE))
write.csv(clas.pop[, c('sim', 'method', 'n.dims', 'fn', 'fp')], 
          file.path("data", "result_summaries", "SIM_Etrue.csv"),
          row.names = FALSE)

# E_{tt, k}
Ett = ddply(clas.tt, 
            c('sim', 'method', 'n.dims'), 
            summarise, 
            fn = sum(!clas.assigned & clas.true),
            fp = sum(clas.assigned & !clas.true))
write.csv(Ett, 
          file.path("data", "result_summaries", "SIM_Ett.csv"),
          row.names = FALSE)

# E_{tt, j, k} (Pi1)
EttPi1 = ddply(loo.alt2.tt, 
               c('sim', 'method', 'n.dims', 'left.out.obs'), 
               summarise,
               fn = sum(!clas.assigned & clas.true),
               fp = sum(clas.assigned & !clas.true))
write.csv(EttPi1, 
          file.path("data", "result_summaries", "SIM_EttPi1.csv"),
          row.names = FALSE)

# E_{tt, j, k} (Pi2)
EttPi2 = ddply(loo.alt1.tt, 
               c('sim', 'method', 'n.dims', 'left.out.obs'), 
               summarise,
               fn = sum(!clas.assigned & clas.true),
               fp = sum(clas.assigned & !clas.true))
write.csv(EttPi2, 
          file.path("data", "result_summaries", "SIM_EttPi2.csv"),
          row.names = FALSE)


# E_{cv, j, k} (Pi1)
loo.alt2.cv = transform(loo.alt2.cv, 
                        fn = as.numeric(!clas.assigned & clas.true),
                        fp = as.numeric(clas.assigned & !clas.true))
write.csv(loo.alt2.cv[, c('sim', 'method', 'n.dims', 'left.out.obs', 'fn', 'fp')], 
          file.path("data", "result_summaries", "SIM_EcvPi1.csv"),
          row.names = FALSE)

# E_{cv, j, k} (Pi2)
loo.alt1.cv = transform(loo.alt1.cv, 
                        fn = as.numeric(!clas.assigned & clas.true),
                        fp = as.numeric(clas.assigned & !clas.true))
write.csv(loo.alt1.cv[, c('sim', 'method', 'n.dims', 'left.out.obs', 'fn', 'fp')], 
          file.path("data", "result_summaries", "SIM_EcvPi2.csv"),
          row.names = FALSE)


# E_{cv, j, i, k} (Pi2)
loo.cv6 = transform(loo.cv6, 
                    fn = as.numeric(!clas.assigned & clas.true),
                    fp = as.numeric(clas.assigned & !clas.true))
write.csv(loo.cv6[, c('sim', 'method', 'n.dims', 'left.out.obs.outer', 'left.out.obs.inner', 'fn', 'fp')], 
          file.path("data", "result_summaries", "SIM_Ecv6Pi1.csv"),
          row.names = FALSE)










# Clear variables
rm(list=ls())

# E_{true, k}
Etrue = read.csv(file.path("output", "result_summaries", "SIM_Etrue.csv"))
Etrue$err = Etrue$fn + Etrue$fp
Etrue[Etrue$err < 1e-16, 'err'] = 0
Etrue = Etrue[, c('sim', 'method', 'n.dims', 'err')]

# E_{tt, k}
Ett = read.csv(file.path("output", "result_summaries", "SIM_Ett.csv"))

# E_{tt, j, k} (Pi1)
EttPi1 = read.csv(file.path("output", "result_summaries", "SIM_EttPi1.csv"))

# E_{tt, j, k} (Pi2)
EttPi2 = read.csv(file.path("output", "result_summaries", "SIM_EttPi2.csv"))

# E_{cv, j, k} (Pi1)
EcvPi1 = read.csv(file.path("output", "result_summaries", "SIM_EcvPi1.csv"))

# E_{cv, j, k} (Pi2)
EcvPi2 = read.csv(file.path("output", "result_summaries", "SIM_EcvPi2.csv"))

# E_{cv, j, i, k} (Pi2)
EcvjikPi1 = read.csv(file.path("output", "result_summaries", "SIM_Ecv6Pi1.csv"))
Ecv6Pi1 = ddply(EcvjikPi1, 
                c('sim', 'method', 'n.dims', 'left.out.obs.outer'), 
                summarise,
                fn = sum(fn), 
                fp = sum(fp))
names(Ecv6Pi1)[names(Ecv6Pi1) == 'left.out.obs.outer'] = 'left.out.obs'
rm('EcvjikPi1')

loo.results = data.frame()
tru.results = data.frame()

# CV1, Pi1 (Alt2)
EcvkPi1 = ddply(EcvPi1, 
                c('sim', 'method', 'n.dims'), 
                summarise,
                err = sum(fn) + sum(fp))

kstarPi1 = ddply(EcvkPi1, 
                 c('sim', 'method'), 
                 summarise,
                 n.dims = min(n.dims[err == min(err)]))

loo.results = rbind(loo.results, transform(merge(EcvkPi1, kstarPi1), n.dims.min = NA, n.dims.max = NA, CV = 1, Pi = 1, Alt = 2))
tru.results = rbind(tru.results, transform(merge(Etrue, kstarPi1), CV = 1, Pi = 1, Alt = 2))


# CV6, Pi1 (Extra Alt)
Ecv6Pi1$err = Ecv6Pi1$fn + Ecv6Pi1$fp
kstarcv6Pi1 = ddply(Ecv6Pi1, 
                    c('sim', 'method', 'left.out.obs'), 
                    summarise,
                    n.dims = min(n.dims[err == min(err)]))

loo.results = rbind(loo.results, transform(ddply(merge(EcvPi1, kstarcv6Pi1), 
                                                 c('sim', 'method'), 
                                                 summarise, 
                                                 err = sum(fn) + sum(fp), 
                                                 n.dims.min = min(n.dims), 
                                                 n.dims.max = max(n.dims)), 
                                           n.dims = NA, CV = 3, Pi = 1, Alt = 5))
tru.results = rbind(tru.results, transform(merge(Etrue, kstarPi1), CV = 3, Pi = 1, Alt = 5))







# CV1, Pi2 (Alt 3)
EcvkPi2 = ddply(EcvPi2, 
                c('sim', 'method', 'n.dims'), 
                summarise,
                err = sum(fn) + sum(fp))

kstarPi2 = ddply(EcvkPi2, 
                 c('sim', 'method'), 
                 summarise,
                 n.dims = min(n.dims[err == min(err)]))

loo.results = rbind(loo.results, transform(merge(EcvkPi2, kstarPi2), n.dims.min = NA, n.dims.max = NA, CV = 1, Pi = 2, Alt = 3))
tru.results = rbind(tru.results, transform(merge(Etrue, kstarPi2), CV = 1, Pi = 2, Alt = 3))


# CV2 kstar_true

# NOTE: min(fn + fp) is zero in every case.
# NOTE: Etrue is identical for Alt 0 and Alt 1
kstar.tru = ddply(Ett, 
                  c('sim', 'method'), 
                  summarise, 
                  n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))

# CV2, Pi1 (Alt 0)
kstar = ddply(EttPi1, 
              c('sim', 'method', 'left.out.obs'), 
              summarise,
              n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results, 
                    ddply(merge(EcvPi1, kstar), 
                          c('sim', 'method'),
                          summarise,
                          err = sum(fn + fp), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 2, 
                          Pi = 1, 
                          Alt = 0))
tru.results = rbind(tru.results, transform(merge(Etrue, kstar.tru), CV = 2, Pi = 1, Alt = 0))

# CV2, Pi2 (Alt 1)
kstar = ddply(EttPi2, 
              c('sim', 'method', 'left.out.obs'), 
              summarise,
              n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results, 
                    ddply(merge(EcvPi2, kstar), 
                          c('sim', 'method'),
                          summarise,
                          err = sum(fn + fp), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 2, 
                          Pi = 2, 
                          Alt = 1))
tru.results = rbind(tru.results, transform(merge(Etrue, kstar.tru), CV = 2, Pi = 2, Alt = 1))
  



# CV3, Pi1 (Alt 4)
kstar = ddply(EcvPi1,
              c('sim', 'method', 'left.out.obs'),
              summarise,
              n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results,
                    ddply(merge(EcvPi1, kstar),
                          c('sim', 'method'),
                          summarise,
                          err = sum(fn + fp),
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 3,
                          Pi = 1,
                          Alt = 4))
tru.results = rbind(tru.results, transform(merge(Etrue, kstarPi1), CV = 3, Pi = 1, Alt = 4))

kstar = ddply(EcvPi2,
              c('sim', 'method', 'left.out.obs'),
              summarise,
              n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results,
                    ddply(merge(EcvPi2, kstar),
                          c('sim', 'method'),
                          summarise,
                          err = sum(fn + fp),
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 3,
                          Pi = 2,
                          Alt = 5))
tru.results = rbind(tru.results, transform(merge(Etrue, kstarPi2), CV = 3, Pi = 2, Alt = 5))








write.csv(loo.results, 
          file.path("data", "result_summaries", "SIM_loo_results.csv"),
          row.names = FALSE)

write.csv(tru.results, 
          file.path("data", "result_summaries", "SIM_tru_results.csv"),
          row.names = FALSE)





















