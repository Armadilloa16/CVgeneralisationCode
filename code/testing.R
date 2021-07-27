
source('housekeeping_functions.R')

# Load clinical variables.
clin.df = load_clin()

# Load patient averages as produced by sumamrise_patients.R
pMeans = load_patient_summary()

# Reshape data from a sparse (long) format into a data-matrix (wide) format.
data.w = reshape_data_matrix(pMeans, clin.df)

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

clas.tt.meta = data.frame()
for (rep in 1:100) {
  clas.tt = data.frame()
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
  names(clas.tt)[names(clas.tt) == 'clas.assigned'] = paste('clas.assigned', rep, sep = '.')
  if (rep == 1) {
    clas.tt.meta = clas.tt
  } else {
    clas.tt.meta = merge(clas.tt.meta, clas.tt)
  }
}