
# Load clinical variables.
clin.df = load_clin()

# Load patient averages as produced by sumamrise_patients.R
pMeans = load_patient_summary()

# Reshape data from a sparse (long) format into a data-matrix (wide) format.
data.w = reshape_data_matrix(pMeans, clin.df)

clas.tt = data.frame()
Ecv.Pi2 = data.frame()
Ett.Pi2 = data.frame()
Ecv.Pi1 = data.frame()
Ett.Pi1 = data.frame()

Ecv.CV3 = data.frame()

# Sample size
n = get_sample_size(clin.df)

# y.var is the class variable
y.var = get_clin_variable_name()

# The class labels.
y = clin.df[match(data.w$Patient, clin.df$Patient), y.var]
# The matrix form for the class labels.
y.m = matrix(FALSE, nrow = n, ncol = 2)
y.m[, 1] = !y
y.m[, 2] = y
# The centred class labels.
y.c = y - mean(y)


# In this case, number of variables ranked by CCA to use.
n.dims = 2:(n-4)



#############################################
### Testing and Training on the full data ###
#############################################
cat("      TT\n")


x.data = data.w[, 2:dim(data.w)[2]]

# Calculate the center (mean) vector of the data not including the
# left-out patient
x.cent = colMeans(x.data)
# Center the data
x.c = x.data - t(matrix(rep(x.cent, n), nrow = dim(x.data)[2]))
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
  
  m = train_lda(x.data[, cca_rank[1:k]], y)
  d = m$d
  cutoff = m$cutoff
  
  tmp = as.matrix(x.data[, cca_rank[1:k]]) %*% d > cutoff
  clas.tt = rbind(clas.tt, data.frame(obs = 1:n,
                                      n.dims = k,
                                      clas.assigned = tmp,
                                      clas.true = y))
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
    
    m = train_lda(x.data[1:n != i, cca_rank[1:k]], y[1:n != i])
    d = m$d
    cutoff = m$cutoff
    
    # CCA
    # Test on the same (n - 1) observations used to train the rule.
    tmp = as.matrix(x.data[1:n != i, cca_rank[1:k]]) %*% d > cutoff
    Ett.Pi2 = rbind(Ett.Pi2, data.frame(left.out.obs = i,
                                        obs = (1:n)[1:n != i],
                                        n.dims = k,
                                        clas.assigned = tmp,
                                        clas.true = y[1:n != i]))
    
    # Apply that rule to the left-out test observation.
    tmp = as.matrix(x.data[1:n == i, cca_rank[1:k]]) %*% d > cutoff
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
  # if (i %% 10 == 0) {
  cat(paste("", i, "\n"))
  # }
  
  # Remove that observation
  x = x.data[1:n != i, ]
  
  # Calculate the center (mean) vector of the data not including the
  # left-out patient
  x.cent = colMeans(x)
  # Center the data
  x.c = x - t(matrix(rep(x.cent, n - 1), nrow = dim(x)[2]))
  # Caluclate the singular value decompostion of the centred data
  x.svd = svd(x.c)
  
  # ################################################################# #
  #   CCA1 - As in Winderbaum et. al. (2016),                         #
  #          essentially phi_1 as in Section 13.3.1 of Koch (2013)    #
  # ################################################################# #
  # Calculate the ranking vector for cca1
  phi1 = x.svd$v %*% diag(x.svd$d^-1) %*% t(x.svd$u) %*% y.c[1:n != i]
  # and the ranks
  cca_rank = order(abs(phi1), decreasing = TRUE)
  
  # Then for each number of variables in n.dims,
  for (k in n.dims) {
    
    m = train_lda(x.data[1:n != i, cca_rank[1:k]], y[1:n != i])
    d = m$d
    cutoff = m$cutoff
    
    # CCA1
    # Test on the same (n - 1) observations used to train the rule.
    tmp = as.matrix(x.data[1:n != i, cca_rank[1:k]]) %*% d > cutoff
    Ett.Pi1 = rbind(Ett.Pi1, data.frame(left.out.obs = i,
                                        obs = (1:n)[1:n != i],
                                        n.dims = k,
                                        clas.assigned = tmp,
                                        clas.true = y[1:n != i]))
    
    # Apply that rule to the left-out test observation.
    tmp = as.matrix(x.data[1:n == i, cca_rank[1:k]]) %*% d > cutoff
    Ecv.Pi1 = rbind(Ecv.Pi1, data.frame(left.out.obs = i,
                                        n.dims = k,
                                        clas.assigned = tmp,
                                        clas.true = y[i]))
  }
  
  # Double CV, i.e. [CV3] 
  for (j in which(1:n != i)) {
    cat(paste("", j))
    
    # Remove the observations
    x = x.data[!(1:n %in% c(i, j)), ]
    
    # Calculate the center (mean) vector of the data not including the
    # left-out patient
    x.cent = colMeans(x)
    # Center the data
    x.c = x - t(matrix(rep(x.cent, n - 1), nrow = dim(x)[2]))
    # Caluclate the singular value decompostion of the centred data
    x.svd = svd(x.c)
    
    # ################################################################# #
    #   CCA1 - As in Winderbaum et. al. (2016),                         #
    #          essentially phi_1 as in Section 13.3.1 of Koch (2013)    #
    # ################################################################# #
    # Calculate the ranking vector for cca1
    phi1 = x.svd$v %*% diag(x.svd$d^-1) %*% t(x.svd$u) %*% y.c[!(1:n %in% c(i, j))]
    # and the ranks
    cca_rank = order(abs(phi1), decreasing = TRUE)
    
    # Then for each number of variables in n.dims,
    for (k in n.dims) {
      
      m = train_lda(x.data[!(1:n %in% c(i, j)), cca_rank[1:k]], y[!(1:n %in% c(i, j))], verbose = TRUE)
      d = m$d
      cutoff = m$cutoff  
      
      # Apply that rule to the left-out test observation.
      tmp = as.matrix(x.data[1:n == j, cca_rank[1:k]]) %*% d > cutoff
      Ecv.CV3 = rbind(Ecv.CV3, data.frame(left.out.obs.outer = i,
                                          left.out.obs.inner = j,
                                          n.dims = k,
                                          clas.assigned = tmp,
                                          clas.true = y[j]))
      
      
    }
  }
}


cat('  Writing output\n')
write.csv(clas.tt, file.path("data", "results", "cca_lda_tt.csv"),
          row.names = FALSE)
write.csv(Ecv.Pi2, file.path("data", "results", "cca_lda_Pi2_Ecv.csv"),
          row.names = FALSE)
write.csv(Ett.Pi2, file.path("data", "results", "cca_lda_Pi2_Ett.csv"),
          row.names = FALSE)
write.csv(Ecv.Pi1, file.path("data", "results", "cca_lda_Pi1_Ecv.csv"),
          row.names = FALSE)
write.csv(Ett.Pi1, file.path("data", "results", "cca_lda_Pi1_Ett.csv"),
          row.names = FALSE)

write.csv(Ecv.CV3, file.path("data", "results", "cca_lda_CV3_Ecv.csv"),
          row.names = FALSE)





