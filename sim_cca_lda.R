
# Load population parameters
mean.p = load_sim_param(TRUE, 'mean')
cov.p  = load_sim_param(TRUE, 'cov')
mean.n = load_sim_param(FALSE, 'mean')
cov.n  = load_sim_param(FALSE, 'cov')

# For each simulation
N = 200
for (sim in 1:N) {
  cat(paste('Simulation:', toString(sim), '\n'))
  
  clas.tt = data.frame()
  clas.pop = data.frame()
  Ecv.Pi2 = data.frame()
  Ett.Pi2 = data.frame()
  Ecv.Pi1 = data.frame()
  Ett.Pi1 = data.frame()
  
  Ecv.CV3 = data.frame()
  
  
  # Load simulated data
  sim.data = load_sim_data(TRUE, sim)
  n.p = nrow(sim.data)
  sim.data = rbind(sim.data, load_sim_data(FALSE, sim))
  n = nrow(sim.data)
  y = rep(c(TRUE, FALSE), c(n.p, n - n.p))
  
  # In this case, number of variables ranked by CCA to use.
  n.dims = 2:(n-4)
  

  
  # The matrix form for the class labels.
  y.m = matrix(FALSE, nrow = n, ncol = 2)
  y.m[, 1] = !y
  y.m[, 2] = y
  # The centred class labels.
  y.c = y - mean(y)
  
  #############################################
  ### Testing and Training on the full data ###
  #############################################
  cat("      TT\n")
  
  # Calculate the center (mean) vector of the data not including the 
  # left-out patient
  x.cent = colMeans(sim.data)
  # Center the data
  x.c = sim.data - t(matrix(rep(x.cent, n), nrow = dim(sim.data)[2]))
  # Caluclate the singular value decompostion of the centred data
  x.svd = svd(x.c)
  
  # ################################################################# #
  #   CCA - As in Winderbaum et. al. (2016),                          # 
  #          essentially phi_1 as in Section 13.3.1 of Koch (2013)    #
  # ################################################################# #
  # Calculate the ranking vector for cca1
  phi1 = x.svd$v %*% diag(x.svd$d^-1) %*% t(x.svd$u) %*% y.c
  # and the ranks
  cca_rank = order(abs(phi1), decreasing = TRUE)

  for (k in n.dims) {

    m = lda_lw(sim.data[, cca_rank[1:k]], y)
    d = m$d
    cutoff = m$cutoff

    tmp = as.matrix(sim.data[, cca_rank[1:k]]) %*% d > cutoff
    clas.tt = rbind(clas.tt, data.frame(obs = 1:n,
                                        n.dims = k,
                                        clas.assigned = tmp,
                                        clas.true = y))

    clas.pop = rbind(clas.pop, data.frame(n.dims = k,
                                          mu.p = as.numeric(mean.p[cca_rank[1:k]] %*% d),
                                          sd.p = sqrt(as.numeric(t(d) %*% cov.p[cca_rank[1:k], cca_rank[1:k]] %*% d)),
                                          mu.n = as.numeric(mean.n[cca_rank[1:k]] %*% d),
                                          sd.n = sqrt(as.numeric(t(d) %*% cov.n[cca_rank[1:k], cca_rank[1:k]] %*% d)),
                                          cutoff = cutoff,
                                          k = k))


  }

  ###########################################################
  ### Pi2 (dimension reduction on full data, then subset) ###
  ###########################################################
  cat("      [Pi2], LO:")

  # For each observation:
  for (i in 1:n) {
    if (i %% 10 == 0) {
      cat(paste("", i))
    }

    # Then for each number of variables
    for (k in n.dims) {

      m = lda_lw(sim.data[1:n != i, cca_rank[1:k]], y[1:n != i])
      d = m$d
      cutoff = m$cutoff

      # CCA1
      # Test on the same (n - 1) observations used to train the rule.
      tmp = as.matrix(sim.data[1:n != i, cca_rank[1:k]]) %*% d > cutoff
      Ett.Pi2 = rbind(Ett.Pi2, data.frame(left.out.obs = i,
                                                    obs = (1:n)[1:n != i],
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[1:n != i]))

      # Apply that rule to the left-out test observation.
      tmp = as.matrix(sim.data[1:n == i, cca_rank[1:k]]) %*% d > cutoff
      Ecv.Pi2 = rbind(Ecv.Pi2, data.frame(left.out.obs = i,
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[i]))
    }
  }
  cat('\n')



  ##############################################
  ### Pi1 (dimension reduction after subset) ###
  ##############################################
  cat("      [Pi1], LO:")

  # For each observation:
  for (i in 1:n) {
    if (i %% 10 == 0) {
      cat(paste("", i))
    }

    # Remove that observation
    x = sim.data[1:n != i, ]

    # Calculate the center (mean) vector of the data not including the
    # left-out patient
    x.cent = colMeans(x)
    # Center the data
    x.c = x - t(matrix(rep(x.cent, n - 1), nrow = dim(x)[2]))
    # Caluclate the singular value decompostion of the centred data
    x.svd = svd(x.c)

    # ################################################################# #
    #   CCA - As in Winderbaum et. al. (2016),                          #
    #          essentially phi_1 as in Section 13.3.1 of Koch (2013)    #
    # ################################################################# #
    # Calculate the ranking vector for cca1
    phi1 = x.svd$v %*% diag(x.svd$d^-1) %*% t(x.svd$u) %*% y.c[1:n != i]
    # and the ranks
    cca_rank = order(abs(phi1), decreasing = TRUE)

    # Then for each number of variables in n.dims,
    for (k in n.dims) {

      m = lda_lw(sim.data[1:n != i, cca_rank[1:k]], y[1:n != i])
      d = m$d
      cutoff = m$cutoff

      # CCA1
      # Test on the same (n - 1) observations used to train the rule.
      tmp = as.matrix(sim.data[1:n != i, cca_rank[1:k]]) %*% d > cutoff
      Ett.Pi1 = rbind(Ett.Pi1, data.frame(left.out.obs = i,
                                                    obs = (1:n)[1:n != i],
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[1:n != i]))

      # Apply that rule to the left-out test observation.
      tmp = as.matrix(sim.data[1:n == i, cca_rank[1:k]]) %*% d > cutoff
      Ecv.Pi1 = rbind(Ecv.Pi1, data.frame(left.out.obs = i,
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[i]))
    }
    
    # Double CV, i.e. [CV3] 
    for (j in which(1:n != i)) {
      # Remove that observation
      x = sim.data[!(1:n %in% c(i, j)), ]
      
      # Calculate the center (mean) vector of the data not including the
      # left-out patient
      x.cent = colMeans(x)
      # Center the data
      x.c = x - t(matrix(rep(x.cent, n - 1), nrow = dim(x)[2]))
      # Caluclate the singular value decompostion of the centred data
      x.svd = svd(x.c)
      
      # ################################################################# #
      #   CCA - As in Winderbaum et. al. (2016),                          #
      #          essentially phi_1 as in Section 13.3.1 of Koch (2013)    #
      # ################################################################# #
      # Calculate the ranking vector for cca1
      phi1 = x.svd$v %*% diag(x.svd$d^-1) %*% t(x.svd$u) %*% y.c[!(1:n %in% c(i, j))]
      # and the ranks
      cca_rank = order(abs(phi1), decreasing = TRUE)
      
      # Then for each number of variables in n.dims,
      for (k in n.dims) {
        
        m = lda_lw(sim.data[!(1:n %in% c(i, j)), cca_rank[1:k]], y[!(1:n %in% c(i, j))])
        d = m$d
        cutoff = m$cutoff
        
        # Apply that rule to the left-out test observation.
        tmp = as.matrix(sim.data[1:n == j, cca_rank[1:k]]) %*% d > cutoff
        Ecv.CV3 = rbind(Ecv.CV3, data.frame(left.out.obs.outer = i,
                                                    left.out.obs.inner = j,
                                                    n.dims = k,
                                                    clas.assigned = tmp,
                                                    clas.true = y[j]))
      }
    }
  }
  
  
      
  cat('  Writing output\n')
  write.csv(clas.tt, file.path("data", "sim_results",
                               paste0("cca_lda_tt_", toString(sim), ".csv")),
            row.names = FALSE)

  write.csv(clas.pop, file.path("data", "sim_results",
                                paste0("cca_lda_pop_", toString(sim), ".csv")),
            row.names = FALSE)

  write.csv(Ecv.Pi2, file.path("data", "sim_results",
                                    paste0("cca_lda_Pi2_Ecv_",
                                           toString(sim), ".csv")),
            row.names = FALSE)

  write.csv(Ett.Pi2, file.path("data", "sim_results",
                                    paste0("cca_lda_Pi2_Ett_",
                                           toString(sim), ".csv")),
            row.names = FALSE)
  write.csv(Ecv.Pi1, file.path("data", "sim_results",
                                    paste0("cca_lda_Pi1_Ecv_",
                                           toString(sim), ".csv")),
            row.names = FALSE)

  write.csv(Ett.Pi1, file.path("data", "sim_results",
                                    paste0("cca_lda_Pi1_Ett_",
                                           toString(sim), ".csv")),
            row.names = FALSE)
  
  
  write.csv(Ecv.CV3, file.path("data", "sim_results",
                                    paste0("cca_lda_CV3_Ecv_",
                                           toString(sim), ".csv")),
            row.names = FALSE)
  
}



# 