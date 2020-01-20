
library(plyr)



clas.tt = transform(read.csv(file.path("data", "results", "cca_lda_tt.csv")), 
                    method = 'CCA-LDA')
clas.tt = rbind(clas.tt,
                transform(read.csv(file.path("data", "results",
                                             "pca_lda_tt.csv")),
                          method = 'PCA-LDA'))

loo.alt1.cv = transform(read.csv(file.path("data", "results", "cca_lda_alt1_loo_cv.csv")), 
                        method = 'CCA-LDA')
loo.alt1.cv = rbind(loo.alt1.cv,
                    transform(read.csv(file.path("data", "results",
                                                 "pca_lda_alt1_loo_cv.csv")),
                              method = 'PCA-LDA'))

loo.alt1.tt = transform(read.csv(file.path("data", "results", "cca_lda_alt1_loo_tt.csv")), 
                        method = 'CCA-LDA')
loo.alt1.tt = rbind(loo.alt1.tt,
                    transform(read.csv(file.path("data", "results",
                                                 "pca_lda_alt1_loo_tt.csv")),
                              method = 'PCA-LDA'))

loo.alt2.cv = transform(read.csv(file.path("data", "results", "cca_lda_alt2_loo_cv.csv")), 
                        method = 'CCA-LDA')
loo.alt2.cv = rbind(loo.alt2.cv,
                    transform(read.csv(file.path("data", "results",
                                                 "pca_lda_alt2_loo_cv.csv")),
                              method = 'PCA-LDA'))

loo.alt2.tt = transform(read.csv(file.path("data", "results", "cca_lda_alt2_loo_tt.csv")), 
                        method = 'CCA-LDA')
loo.alt2.tt = rbind(loo.alt2.tt,
                    transform(read.csv(file.path("data", "results",
                                                 "pca_lda_alt2_loo_tt.csv")),
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
EttPi1 = ddply(loo.alt2.tt, 
               c('method', 'n.dims', 'left.out.obs'), 
               summarise,
               fn = sum(!clas.assigned & clas.true),
               fp = sum(clas.assigned & !clas.true))
write.csv(EttPi1, 
          file.path("data", "result_summaries", "ETMA_EttPi1.csv"),
          row.names = FALSE)

# E_{tt, j, k} (Pi2)
EttPi2 = ddply(loo.alt1.tt, 
               c('method', 'n.dims', 'left.out.obs'), 
               summarise,
               fn = sum(!clas.assigned & clas.true),
               fp = sum(clas.assigned & !clas.true))
write.csv(EttPi2, 
          file.path("data", "result_summaries", "ETMA_EttPi2.csv"),
          row.names = FALSE)


# E_{cv, j, k} (Pi1)
loo.alt2.cv = transform(loo.alt2.cv, 
                        fn = as.numeric(!clas.assigned & clas.true),
                        fp = as.numeric(clas.assigned & !clas.true))
EcvPi1 = loo.alt2.cv[, c('method', 'n.dims', 'left.out.obs', 'fn', 'fp')]
write.csv(EcvPi1, 
          file.path("data", "result_summaries", "ETMA_EcvPi1.csv"),
          row.names = FALSE)

# E_{cv, j, k} (Pi2)
loo.alt1.cv = transform(loo.alt1.cv, 
                        fn = as.numeric(!clas.assigned & clas.true),
                        fp = as.numeric(clas.assigned & !clas.true))
EcvPi2 = loo.alt1.cv[, c('method', 'n.dims', 'left.out.obs', 'fn', 'fp')]
write.csv(EcvPi2, 
          file.path("data", "result_summaries", "ETMA_EcvPi2.csv"),
          row.names = FALSE)




loo.results = data.frame()

# CV1, Pi1 (Alt2)
Ecvk = ddply(EcvPi1, 
             c('method', 'n.dims'), 
             summarise,
             err = sum(fn) + sum(fp))
             
kstar = ddply(Ecvk, 
              'method', 
              summarise,
              n.dims = min(n.dims[err == min(err)]))

loo.results = rbind(loo.results, transform(merge(Ecvk, kstar), n.dims.min = NA, n.dims.max = NA, CV = 1, Pi = 1, Alt = 2))

# CV1, Pi2 (Alt 3)
Ecvk = ddply(loo.alt1.cv, 
             c('method', 'n.dims'), 
             summarise,
             err = sum(fn) + sum(fp))

kstar = ddply(Ecvk, 
              'method', 
              summarise,
              n.dims = min(n.dims[err == min(err)]))

loo.results = rbind(loo.results, transform(merge(Ecvk, kstar), n.dims.min = NA, n.dims.max = NA, CV = 1, Pi = 2, Alt = 3))

# CV2, Pi1 (Alt 0)
kstar = ddply(EttPi1, 
              c('method', 'left.out.obs'), 
              summarise,
              n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results, 
                    ddply(merge(EcvPi1, kstar), 
                          'method',
                          summarise,
                          err = sum(fn + fp), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 2, 
                          Pi = 1, 
                          Alt = 0))
                          
# CV2, Pi2 (Alt 1)
kstar = ddply(EttPi2, 
              c('method', 'left.out.obs'), 
              summarise,
              n.dims = min(n.dims[(fn + fp) == min(fn + fp)]))
loo.results = rbind(loo.results, 
                    ddply(merge(EcvPi2, kstar), 
                          'method',
                          summarise,
                          err = sum(fn + fp), 
                          n.dims.min = min(n.dims),
                          n.dims.max = max(n.dims),
                          n.dims = NA,
                          CV = 2, 
                          Pi = 2, 
                          Alt = 1))

write.csv(loo.results, 
          file.path("data", "result_summaries", "ETMA_loo_results.csv"),
          row.names = FALSE)

# Note that classification results can be obtained as follows: 
# kstar = ddply(Ett, 'method', summarise, n.dims = min(n.dims[(fp + fn) == min(fp + fn)]))
# merge(kstar, Ett)



