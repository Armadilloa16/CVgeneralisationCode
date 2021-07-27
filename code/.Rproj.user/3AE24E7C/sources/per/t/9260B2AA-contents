
library(plyr)

clas.tt = transform(read.csv(file.path("data", "results", "cca_lda_tt.csv")), 
                    method = 'CCA-LDA')
clas.tt = rbind(clas.tt,
                transform(read.csv(file.path("data", "results", 
                                             "pca_lda_tt.csv")),
                          method = 'PCA-LDA'))

Ecv.Pi2 = transform(read.csv(file.path("data", "results", 
                                       "cca_lda_Pi2_Ecv.csv")), 
                        method = 'CCA-LDA')
Ecv.Pi2 = rbind(Ecv.Pi2,
                    transform(read.csv(file.path("data", "results",
                                                 "pca_lda_Pi2_Ecv.csv")),
                              method = 'PCA-LDA'))

Ett.Pi2 = transform(read.csv(file.path("data", "results", 
                                       "cca_lda_Pi2_Ett.csv")), 
                        method = 'CCA-LDA')
Ett.Pi2 = rbind(Ett.Pi2,
                    transform(read.csv(file.path("data", "results",
                                                 "pca_lda_Pi2_Ett.csv")),
                              method = 'PCA-LDA'))

Ecv.Pi1 = transform(read.csv(file.path("data", "results", 
                                       "cca_lda_Pi1_Ecv.csv")), 
                        method = 'CCA-LDA')
Ecv.Pi1 = rbind(Ecv.Pi1,
                    transform(read.csv(file.path("data", "results",
                                                 "pca_lda_Pi1_Ecv.csv")),
                              method = 'PCA-LDA'))

Ett.Pi1 = transform(read.csv(file.path("data", "results", 
                                       "cca_lda_Pi1_Ett.csv")), 
                        method = 'CCA-LDA')
Ett.Pi1 = rbind(Ett.Pi1,
                    transform(read.csv(file.path("data", "results",
                                                 "pca_lda_Pi1_Ett.csv")),
                              method = 'PCA-LDA'))


Ecv.CV3 = transform(read.csv(file.path("data", "results", 
                                       "cca_lda_CV3_Ecv.csv")), 
                    method = 'CCA-LDA')
Ecv.CV3 = rbind(Ecv.CV3,
                transform(read.csv(file.path("data", "results",
                                             "pca_lda_CV3_Ecv.csv")),
                          method = 'PCA-LDA'))

# E_{tt, k}
Ett = ddply(clas.tt, 
            c('method', 'n.dims'), 
            summarise, 
            fn = sum(!clas.assigned & clas.true),
            fp = sum(clas.assigned & !clas.true))
write.csv(Ett, 
          file.path("data", "result_summaries", "ETMA_Ett.csv"),
          row.names = FALSE)

# E_{tt, j, k} (Pi1)
Ettj.Pi1 = ddply(Ett.Pi1, 
               c('method', 'n.dims', 'left.out.obs'), 
               summarise,
               fn = sum(!clas.assigned & clas.true),
               fp = sum(clas.assigned & !clas.true))
write.csv(Ettj.Pi1, 
          file.path("data", "result_summaries", "ETMA_Ettj_Pi1.csv"),
          row.names = FALSE)

# E_{tt, j, k} (Pi2)
Ettj.Pi2 = ddply(Ett.Pi2, 
               c('method', 'n.dims', 'left.out.obs'), 
               summarise,
               fn = sum(!clas.assigned & clas.true),
               fp = sum(clas.assigned & !clas.true))
write.csv(Ettj.Pi2, 
          file.path("data", "result_summaries", "ETMA_Ettj_Pi2.csv"),
          row.names = FALSE)


# E_{cv, j, k} (Pi1)
Ecvj.Pi1 = transform(Ecv.Pi1, 
                        fn = as.numeric(!clas.assigned & clas.true),
                        fp = as.numeric(clas.assigned & !clas.true))
Ecvj.Pi1 = Ecvj.Pi1[, c('method', 'n.dims', 'left.out.obs', 'fn', 'fp')]
write.csv(Ecvj.Pi1, 
          file.path("data", "result_summaries", "ETMA_Ecvj_Pi1.csv"),
          row.names = FALSE)

# E_{cv, j, k} (Pi2)
Ecvj.Pi2 = transform(Ecv.Pi2, 
                        fn = as.numeric(!clas.assigned & clas.true),
                        fp = as.numeric(clas.assigned & !clas.true))
Ecvj.Pi2 = Ecvj.Pi2[, c('method', 'n.dims', 'left.out.obs', 'fn', 'fp')]
write.csv(Ecvj.Pi2, 
          file.path("data", "result_summaries", "ETMA_Ecvj_Pi2.csv"),
          row.names = FALSE)

# E_{cv, j, k} (CV3)
Ecvj.CV3 = ddply(Ecv.CV3, 
                 c('method', 'n.dims', 'left.out.obs.outer'), 
                 summarise, 
                 fn = sum(as.numeric(!clas.assigned & clas.true)),
                 fp = sum(as.numeric(clas.assigned & !clas.true)))
names(Ecvj.CV3)[names(Ecvj.CV3) == 'left.out.obs.outer'] = 'left.out.obs'
write.csv(Ecvj.CV3, 
          file.path("data", "result_summaries", "ETMA_Ecvj_CV3.csv"),
          row.names = FALSE)





rm(list = ls())
Ett = read.csv(file.path("data", "result_summaries", "ETMA_Ett.csv"))
Ettj.Pi1 = read.csv(file.path("data", "result_summaries", "ETMA_Ettj_Pi1.csv"))
Ettj.Pi2 = read.csv(file.path("data", "result_summaries", "ETMA_Ettj_Pi2.csv"))
Ecvj.Pi1 = read.csv(file.path("data", "result_summaries", "ETMA_Ecvj_Pi1.csv"))
Ecvj.Pi2 = read.csv(file.path("data", "result_summaries", "ETMA_Ecvj_Pi2.csv"))
Ecvj.CV3 = read.csv(file.path("data", "result_summaries", "ETMA_Ecvj_CV3.csv"))






# E_{cv, k} (Pi1)
Ecv.Pi1 = ddply(Ecvj.Pi1, 
                c('method', 'n.dims'), 
                summarise,
                err = sum(fn) + sum(fp))

# E_{cv, k} (Pi2)
Ecv.Pi2 = ddply(Ecvj.Pi2, 
                 c('method', 'n.dims'), 
                 summarise,
                 err = sum(fn) + sum(fp))




loo.results = data.frame()

# [CV1], [Pi1]
kstar = ddply(Ecv.Pi1, 
              'method', 
              summarise,
              n.dims = min(n.dims[err == min(err)]))
loo.results = rbind(loo.results, transform(merge(Ecv.Pi1, kstar), 
                                           n.dims.min = NA, n.dims.max = NA, 
                                           CV = 1, Pi = 1))

# [CV1], [Pi2]
kstar = ddply(Ecv.Pi2, 
              'method', 
              summarise,
              n.dims = min(n.dims[err == min(err)]))
loo.results = rbind(loo.results, transform(merge(Ecv.Pi2, kstar), 
                                           n.dims.min = NA, n.dims.max = NA, 
                                           CV = 1, Pi = 2))

# [CV2], [Pi1]
kstarj = ddply(Ettj.Pi1, 
               c('method', 'left.out.obs'), 
               summarise,
               n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results, 
                    ddply(merge(Ecvj.Pi1, kstarj), 
                          'method',
                          summarise,
                          err = sum(fn + fp), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 2, 
                          Pi = 1))
                          
# [CV2], [Pi2]
kstarj = ddply(Ettj.Pi2, 
               c('method', 'left.out.obs'), 
               summarise,
               n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results, 
                    ddply(merge(Ecvj.Pi2, kstarj), 
                          'method',
                          summarise,
                          err = sum(fn + fp), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 2, 
                          Pi = 2))

# [CV3]
kstarj = ddply(Ecvj.CV3, 
              c('method', 'left.out.obs'), 
              summarise,
              n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results, 
                    ddply(merge(Ecvj.Pi1, kstarj), 
                          'method',
                          summarise,
                          err = sum(fn + fp), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 3, 
                          Pi = 1))


write.csv(loo.results, 
          file.path("data", "result_summaries", "ETMA_loo_results.csv"),
          row.names = FALSE)

# Note that classification results can be obtained as follows: 
# kstar = ddply(Ett, 'method', summarise, n.dims = min(n.dims[(fp + fn) == min(fp + fn)]))
# merge(kstar, Ett)



