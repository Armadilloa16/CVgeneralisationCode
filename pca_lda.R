
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
# The centred class labels.
y.c = y - mean(y)

x.data = data.w[, 2:dim(data.w)[2]]

# Note that PCA scores can be reproduced by centering and multipling by the 
# eigenvectors (rotation matrix) in the following way, if we let:
x = prcomp(x.data)$x

# In this case, number of PCs to use.
n.dims = 2:(n-4)

#############################################
### Testing and Training on the full data ###
#############################################
cat("      TT\n")

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
  
  # Then for each number of principal components in n.dims,
  for (k in n.dims) {
    # Fit a Fishers Linear Discriminant Analysis model to the Training data
    # using the first k principal component scores as the variables.
    m = train_lda(x[1:n != i, 1:k], y[1:n != i])
    d = m$d
    cutoff = m$cutoff
    
    # Test on the same (n - 1) observations used to train the rule.
    tmp = x[1:n != i, 1:k] %*% d > cutoff
    Ett.Pi2 = rbind(Ett.Pi2, data.frame(left.out.obs = i,
                                        obs = (1:n)[1:n != i],
                                        n.dims = k,
                                        clas.assigned = tmp,
                                        clas.true = y[1:n != i]))
    
    # Apply that rule to the left-out test observation.
    tmp = x[1:n == i, 1:k] %*% d > cutoff
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
  x = x.data[1:n != i, ]
  # calculate the principal components of training data
  pca = prcomp(x)
  # and apply the same centering and rotation to the test observation.
  z = (as.matrix(x.data[1:n == i, ]) - pca$center) %*% pca$rotation
  
  # Then for each number of principal components in n.dims,
  for (k in n.dims) {
    # Fit a Fishers Linear Discriminant Analysis model to the Training data
    # using the first k principal component scores as the variables.
    m = train_lda(pca$x[, 1:k], y[1:n != i])
    d = m$d
    cutoff = m$cutoff
    
    # Test on the same (n - 1) observations used to train the rule.
    tmp = pca$x[, 1:k] %*% d > cutoff
    Ett.Pi1 = rbind(Ett.Pi1, data.frame(left.out.obs = i,
                                        obs = (1:n)[1:n != i],
                                        n.dims = k,
                                        clas.assigned = tmp,
                                        clas.true = y[1:n != i]))
    
    # Apply that rule to the left-out test observation z.
    tmp = z[1:k] %*% d > cutoff
    Ecv.Pi1 = rbind(Ecv.Pi1, data.frame(left.out.obs = i,
                                        n.dims = k,
                                        clas.assigned = tmp,
                                        clas.true = y[i]))
  }
  
  # Double CV, i.e. [CV3] 
  for (j in which(1:n != i)) {
    # Remove both observations
    x = x.data[!(1:n %in% c(i, j)), ]
    # calculate the principal components of training data
    pca = prcomp(x)
    # and apply the same centering and rotation to the test observation.
    z = (as.matrix(x.data[1:n == j, ]) - pca$center) %*% pca$rotation
    
    # Then for each number of principal components in n.dims,
    for (k in n.dims) {
      # Fit a Fishers Linear Discriminant Analysis model to the Training data
      # using the first k principal component scores as the variables.
      m = train_lda(pca$x[, 1:k], y[!(1:n %in% c(i, j))])
      d = m$d
      cutoff = m$cutoff
      
      # Apply that rule to the left-out test observation z.
      tmp = z[1:k] %*% d > cutoff
      Ecv.CV3 = rbind(Ecv.CV3, data.frame(left.out.obs.outer = i,
                                          left.out.obs.inner = j,
                                          n.dims = k,
                                          clas.assigned = tmp,
                                          clas.true = y[j]))
      
    }
  }
}



cat('  Writing output\n')
write.csv(clas.tt, file.path("data", "results", "pca_lda_tt.csv"),
          row.names = FALSE)
write.csv(Ecv.Pi2, file.path("data", "results", "pca_lda_Pi2_Ecv.csv"),
          row.names = FALSE)
write.csv(Ett.Pi2, file.path("data", "results", "pca_lda_Pi2_Ett.csv"),
          row.names = FALSE)
write.csv(Ecv.Pi1, file.path("data", "results", "pca_lda_Pi1_Ecv.csv"),
          row.names = FALSE)
write.csv(Ett.Pi1, file.path("data", "results", "pca_lda_Pi1_Ett.csv"),
          row.names = FALSE)

write.csv(Ecv.CV3, file.path("data", "results", "pca_lda_CV3_Ecv.csv"),
          row.names = FALSE)

