
# Load population parameters
mean.p.pop = load_sim_param(TRUE, 'mean')
cov.p.pop  = load_sim_param(TRUE, 'cov')
mean.n.pop = load_sim_param(FALSE, 'mean')
cov.n.pop  = load_sim_param(FALSE, 'cov')

# For each simulation
N = 200
for (sim in 1:N) {
  cat(paste('Simulation:', toString(sim), '\n'))
  
  clas.tt = data.frame()
  clas.pop = data.frame()
  clas.df.alt1 = data.frame()
  misc.df.alt1 = data.frame()
  clas.df.alt2 = data.frame()
  misc.df.alt2 = data.frame()
  
  clas.df.cv6 = data.frame()
  
  # Load simulated data
  sim.data = load_sim_data(TRUE, sim)
  n.p = nrow(sim.data)
  sim.data = rbind(sim.data, load_sim_data(FALSE, sim))
  n = nrow(sim.data)
  y = rep(c(TRUE, FALSE), c(n.p, n - n.p))
  
  # In this case, number of PCs to use.
  n.dims = 2:(n-4)
  
  
  #############################################
  ### Testing and Training on the full data ###
  #############################################
  cat("      TT\n")

  # Note that PCA scores can be reproduced by centering and multipling by the
  # eigenvectors (rotation matrix) in the following way, if we let:
  pca = prcomp(as.matrix(sim.data))
  x = pca$x

  # Project population parameters
  mean.p = (mean.p.pop - pca$center) %*% pca$rotation
  cov.p  = t(pca$rotation) %*% cov.p.pop %*% pca$rotation
  mean.n = (mean.n.pop - pca$center) %*% pca$rotation
  cov.n  = t(pca$rotation) %*% cov.n.pop %*% pca$rotation

  # Then the principal component scores as in pca$x are the same as
  #
  # ((as.matrix(x) - t(matrix(rep(pca$center,n), ncol = n))) %*% pca$rotation)

  # Then for each number of principal components in n.dims,
  for (k in n.dims) {
    # Fit a Fishers Linear Discriminant Analysis model to the Training data
    # using the first k principal component scores as the variables.
    m  = train_lda(x[, 1:k], y)
    d = m$d
    cutoff = m$cutoff
    tmp = x[, 1:k] %*% d > cutoff
    clas.tt = rbind(clas.tt, data.frame(obs = 1:n,
                                        n.dims = k,
                                        clas.assigned = tmp,
                                        clas.true = y))

    clas.pop = rbind(clas.pop, data.frame(n.dims = k,
                                          mu.p = as.numeric(mean.p[1:k] %*% d),
                                          sd.p = sqrt(as.numeric(t(d) %*% cov.p[1:k, 1:k] %*% d)),
                                          mu.n = as.numeric(mean.n[1:k] %*% d),
                                          sd.n = sqrt(as.numeric(t(d) %*% cov.n[1:k, 1:k] %*% d)),
                                          cutoff = cutoff,
                                          k = k))
  }
  
  ###################################
  ### Alternatives 1 and 3 LOO-CV ###
  ###################################
  cat("      Alt1, LO:")

  # For each observation:
  for (i in 1:n) {
    if (i %% 10 == 0) {
      cat(paste("", i))
    }

    # Then for each number of principal components in n.dims,
    for (k in n.dims) {
      # Fit a Fishers Linear Discriminant Analysis model to the Training data
      # using the first k principal component scores as the variables.
      m  = train_lda(x[1:n != i, 1:k], y[1:n != i])
      d = m$d
      cutoff = m$cutoff

      # Test on the same (n - 1) observations used to train the rule.
      tmp = x[1:n != i, 1:k] %*% d > cutoff
      misc.df.alt1 = rbind(misc.df.alt1, data.frame(left.out.obs = i,
                                                    obs = (1:n)[1:n != i],
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[1:n != i]))

      # Apply that rule to the left-out test observation.
      tmp = x[1:n == i, 1:k] %*% d > cutoff
      clas.df.alt1 = rbind(clas.df.alt1, data.frame(left.out.obs = i,
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[i]))
    }
  }
  cat('\n')



  #######################################################
  ### Current Implementation and Alternative 2 LOO-CV ###
  #######################################################
  cat("      Alt2, LO:")

  # For each observation:
  for (i in 1:n) {
    if (i %% 10 == 0) {
      cat(paste("", i))
    }

    # Remove that observation
    x = sim.data[1:n != i, ]
    # calculate the principal components of training data
    pca = prcomp(x)
    # and apply the same centering and rotation to the test observation.
    z = (as.matrix(sim.data[1:n == i, ]) - pca$center) %*% pca$rotation

    # Then for each number of principal components in n.dims,
    for (k in n.dims) {
      # Fit a Fishers Linear Discriminant Analysis model to the Training data
      # using the first k principal component scores as the variables.
      m = train_lda(pca$x[, 1:k], y[1:n != i])
      d = m$d
      cutoff = m$cutoff

      # Test on the same (n - 1) observations used to train the rule.
      tmp = pca$x[, 1:k] %*% d > cutoff
      misc.df.alt2 = rbind(misc.df.alt2, data.frame(left.out.obs = i,
                                                    obs = (1:n)[1:n != i],
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[1:n != i]))

      # Apply that rule to the left-out test observation z.
      tmp = z[1:k] %*% d > cutoff
      clas.df.alt2 = rbind(clas.df.alt2, data.frame(left.out.obs = i,
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[i]))
    }
    
    
    # CV6
    for (j in which(1:n != i)) {
      # Remove the observations
      x = sim.data[!(1:n %in% c(i, j)), ]
      # calculate the principal components of training data
      pca = prcomp(x)
      # and apply the same centering and rotation to the test observation.
      z = (as.matrix(sim.data[1:n == j, ]) - pca$center) %*% pca$rotation
      
      # Then for each number of principal components in n.dims,
      for (k in n.dims) {
        # Fit a Fishers Linear Discriminant Analysis model to the Training data
        # using the first k principal component scores as the variables.
        m = train_lda(pca$x[, 1:k], y[!(1:n %in% c(i, j))])
        d = m$d
        cutoff = m$cutoff
        
        # Apply that rule to the left-out test observation z.
        tmp = z[1:k] %*% d > cutoff
        clas.df.cv6 = rbind(clas.df.cv6, data.frame(left.out.obs.outer = i,
                                                    left.out.obs.inner = j,
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[j]))
        
      }
    }
  }
  
  cat('  Writing output\n')
  write.csv(clas.tt, file.path("data", "sim_results",
                               paste0("pca_lda_tt_", toString(sim), ".csv")),
            row.names = FALSE)

  write.csv(clas.pop, file.path("data", "sim_results",
                                paste0("pca_lda_pop_", toString(sim), ".csv")),
            row.names = FALSE)

  write.csv(clas.df.alt1, file.path("data", "sim_results",
                                    paste0("pca_lda_alt1_loo_cv_",
                                           toString(sim), ".csv")),
            row.names = FALSE)

  write.csv(misc.df.alt1, file.path("data", "sim_results",
                                    paste0("pca_lda_alt1_loo_tt_",
                                           toString(sim), ".csv")),
            row.names = FALSE)
  write.csv(clas.df.alt2, file.path("data", "sim_results",
                                    paste0("pca_lda_alt2_loo_cv_",
                                           toString(sim), ".csv")),
            row.names = FALSE)

  write.csv(misc.df.alt2, file.path("data", "sim_results",
                                    paste0("pca_lda_alt2_loo_tt_",
                                           toString(sim), ".csv")),
            row.names = FALSE)
  
  
  write.csv(clas.df.cv6, file.path("data", "sim_results",
                                    paste0("pca_lda_cv6_loo_cv_",
                                           toString(sim), ".csv")),
            row.names = FALSE)
}



