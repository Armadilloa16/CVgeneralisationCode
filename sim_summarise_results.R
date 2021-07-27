
library(plyr)

N = 200

# For each simlation, load the results.
clas.pop = data.frame()
clas.tt = data.frame()
loo.Pi2.cv = data.frame()
loo.Pi2.tt = data.frame()
loo.Pi1.cv = data.frame()
loo.Pi1.tt = data.frame()
loo.CV3 = data.frame()
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



  loo.Pi2.cv = rbind(loo.Pi2.cv,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("cca_lda_Pi2_Ecv_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim, 
                                method = 'CCA-LDA'))
  loo.Pi2.cv = rbind(loo.Pi2.cv,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("pca_lda_Pi2_Ecv_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim,
                                method = 'PCA-LDA'))



  loo.Pi2.tt = rbind(loo.Pi2.tt,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("cca_lda_Pi2_Ett_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim, 
                                method = 'CCA-LDA'))
  loo.Pi2.tt = rbind(loo.Pi2.tt,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("pca_lda_Pi2_Ett_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim,
                                method = 'PCA-LDA'))



  loo.Pi1.cv = rbind(loo.Pi1.cv,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("cca_lda_Pi1_Ecv_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim, 
                                method = 'CCA-LDA'))
  loo.Pi1.cv = rbind(loo.Pi1.cv,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("pca_lda_Pi1_Ecv_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim,
                                method = 'PCA-LDA'))



  loo.Pi1.tt = rbind(loo.Pi1.tt,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("cca_lda_Pi1_Ett_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim, 
                                method = 'CCA-LDA'))
  loo.Pi1.tt = rbind(loo.Pi1.tt,
                      transform(read.csv(file.path("data", "sim_results",
                                                   paste0("pca_lda_Pi1_Ett_",
                                                          toString(sim),
                                                          ".csv"))),
                                sim = sim,
                                method = 'PCA-LDA'))


  
  
  
  loo.CV3 = rbind(loo.CV3,
                  transform(read.csv(file.path("data", "sim_results",
                                               paste0("cca_lda_CV3_Ecv_",
                                                      toString(sim),
                                                      ".csv"))),
                            sim = sim, 
                            method = 'CCA-LDA'))
  loo.CV3 = rbind(loo.CV3,
                  transform(read.csv(file.path("data", "sim_results",
                                               paste0("pca_lda_CV3_Ecv_",
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
Ettj.Pi1 = ddply(loo.Pi1.tt, 
               c('sim', 'method', 'n.dims', 'left.out.obs'), 
               summarise,
               fn = sum(!clas.assigned & clas.true),
               fp = sum(clas.assigned & !clas.true))
write.csv(Ettj.Pi1, 
          file.path("data", "result_summaries", "SIM_Ettj_Pi1.csv"),
          row.names = FALSE)

# E_{tt, j, k} (Pi2)
Ettj.Pi2 = ddply(loo.Pi2.tt, 
               c('sim', 'method', 'n.dims', 'left.out.obs'), 
               summarise,
               fn = sum(!clas.assigned & clas.true),
               fp = sum(clas.assigned & !clas.true))
write.csv(Ettj.Pi2, 
          file.path("data", "result_summaries", "SIM_Ettj_Pi2.csv"),
          row.names = FALSE)


# E_{cv, j, k} (Pi1)
Ecvj.Pi1 = transform(loo.Pi1.cv, 
                        fn = as.numeric(!clas.assigned & clas.true),
                        fp = as.numeric(clas.assigned & !clas.true))
write.csv(Ecvj.Pi1[, c('sim', 'method', 'n.dims', 'left.out.obs', 'fn', 'fp')], 
          file.path("data", "result_summaries", "SIM_Ecvj_Pi1.csv"),
          row.names = FALSE)

# E_{cv, j, k} (Pi2)
Ecvj.Pi2 = transform(loo.Pi2.cv, 
                        fn = as.numeric(!clas.assigned & clas.true),
                        fp = as.numeric(clas.assigned & !clas.true))
write.csv(Ecvj.Pi2[, c('sim', 'method', 'n.dims', 'left.out.obs', 'fn', 'fp')], 
          file.path("data", "result_summaries", "SIM_Ecvj_Pi2.csv"),
          row.names = FALSE)


# E_{cv, j, k} [CV3] (Pi1) 
Ecvj.CV3 = ddply(loo.CV3, 
                 c('sim', 'method', 'n.dims', 'left.out.obs.outer'), 
                 summarise, 
                 fn = sum(as.numeric(!clas.assigned & clas.true)),
                 fp = sum(as.numeric(clas.assigned & !clas.true)))
names(Ecvj.CV3)[names(Ecvj.CV3) == 'left.out.obs.outer'] = 'left.out.obs'
write.csv(Ecvj.CV3, 
          file.path("data", "result_summaries", "SIM_Ecvj_CV3.csv"),
          row.names = FALSE)










# Clear variables
rm(list=ls())

fnfp2err = function(x, tol = 1e-15){
  x$err = x$fn + x$fp
  x[x$err < tol, 'err'] = 0
  x$fn = NULL
  x$fp = NULL
  return(x)
}

# E_{true, k}
Etrue = fnfp2err(read.csv(file.path("data", "result_summaries", "SIM_Etrue.csv")))

# E_{tt, k}
Ett = fnfp2err(read.csv(file.path("data", "result_summaries", "SIM_Ett.csv")))

# E_{tt, j, k} (Pi1)
Ettj.Pi1 = fnfp2err(read.csv(file.path("data", "result_summaries", "SIM_Ettj_Pi1.csv")))

# E_{tt, j, k} (Pi2)
Ettj.Pi2 = fnfp2err(read.csv(file.path("data", "result_summaries", "SIM_Ettj_Pi2.csv")))

# E_{cv, j, k} (Pi1)
Ecvj.Pi1 = fnfp2err(read.csv(file.path("data", "result_summaries", "SIM_Ecvj_Pi1.csv")))
Ecv.Pi1 = ddply(Ecvj.Pi1, 
                c('sim', 'method', 'n.dims'), 
                summarise,
                err = sum(err))

# E_{cv, j, k} (Pi2)
Ecvj.Pi2 = fnfp2err(read.csv(file.path("data", "result_summaries", "SIM_Ecvj_Pi2.csv")))
Ecv.Pi2 = ddply(Ecvj.Pi2, 
                c('sim', 'method', 'n.dims'), 
                summarise,
                err = sum(err))

# E_{cv, j, i, k} (Pi2)
Ecvj.CV3 = fnfp2err(read.csv(file.path("data", "result_summaries", "SIM_Ecvj_CV3.csv")))





loo.results = data.frame()
tru.results = data.frame()

# CV1, Pi1
kstar = ddply(Ecv.Pi1, 
              c('sim', 'method'), 
              summarise,
              n.dims = min(n.dims[err == min(err)]))
tru.results = rbind(tru.results, transform(merge(Etrue, kstar), CV = 1, Pi = 1))
# Note Prediction or true error is identical for CV3, Pi1.
tru.results = rbind(tru.results, transform(merge(Etrue, kstar), CV = 3, Pi = 1))
loo.results = rbind(loo.results, transform(merge(Ecv.Pi1, kstar), 
                                           n.dims.min = NA, n.dims.max = NA, 
                                           CV = 1, Pi = 1))

# CV1, Pi2
kstar = ddply(Ecv.Pi2, 
              c('sim', 'method'), 
              summarise,
              n.dims = min(n.dims[err == min(err)]))
tru.results = rbind(tru.results, transform(merge(Etrue, kstar), CV = 1, Pi = 2))
loo.results = rbind(loo.results, transform(merge(Ecv.Pi2, kstar), 
                                           n.dims.min = NA, n.dims.max = NA, 
                                           CV = 1, Pi = 2))

# CV2, Pi1
# Prediction (population or true error):
kstar = ddply(Ett, 
              c('sim', 'method'), 
              summarise, 
              n.dims = min(n.dims[err == min(err)]))
tru.results = rbind(tru.results, transform(merge(Etrue, kstar), CV = 2, Pi = 1))
# Note Prediction or true error is identical for CV2, Pi2.
tru.results = rbind(tru.results, transform(merge(Etrue, kstar), CV = 2, Pi = 2))
# Estimation (cross validation error):
kstar = ddply(Ettj.Pi1, 
              c('sim', 'method', 'left.out.obs'), 
              summarise,
              n.dims = min(n.dims[err == min(err)]))
loo.results = rbind(loo.results, 
                    ddply(merge(Ecvj.Pi1, kstar), 
                          c('sim', 'method'),
                          summarise,
                          err = sum(err), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 2, 
                          Pi = 1))

# CV2, Pi2
# Note Prediction or true error is identical for CV2, Pi1.
# Estimation (cross validation error):
kstar = ddply(Ettj.Pi2, 
              c('sim', 'method', 'left.out.obs'), 
              summarise,
              n.dims = min(n.dims[err == min(err)]))
loo.results = rbind(loo.results, 
                    ddply(merge(Ecvj.Pi2, kstar), 
                          c('sim', 'method'),
                          summarise,
                          err = sum(err), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 2, 
                          Pi = 2))

# CV3 (double CV)
# Note Prediction or true error is identical for CV1, Pi1.
# Estimation (cross validation error):
kstar = ddply(Ecvj.CV3, 
              c('sim', 'method', 'left.out.obs'), 
              summarise,
              n.dims = min(n.dims[err == min(err)]))
loo.results = rbind(loo.results, 
                    ddply(merge(Ecvj.Pi1, kstar), 
                          c('sim', 'method'), 
                          summarise, 
                          err = sum(err), 
                          n.dims.min = min(n.dims), 
                          n.dims.max = max(n.dims), 
                          n.dims = NA, 
                          CV = 3, 
                          Pi = 1))









write.csv(loo.results, 
          file.path("data", "result_summaries", "SIM_loo_results.csv"),
          row.names = FALSE)

write.csv(tru.results, 
          file.path("data", "result_summaries", "SIM_tru_results.csv"),
          row.names = FALSE)





















