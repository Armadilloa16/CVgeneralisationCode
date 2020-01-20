
# Load clinical variables.
clin.df = load_clin()

# Load patient averages as produced by sumamrise_patients.R
pMeans = load_patient_summary()

# Reshape data from a sparse (long) format into a data-matrix (wide) format.
data.w = reshape_data_matrix(pMeans, clin.df)

# y.var is the class variable, which varies depending on experiement.
y.var = clin_variable_name()

# The class labels.
y = clin.df[match(data.w$Patient, clin.df$Patient), y.var]
sam = data.frame(p = sum(y), n = sum(!y))

mean.p = colMeans(data.w[y,  2:dim(data.w)[2]])
cov.p  =      cov(data.w[y,  2:dim(data.w)[2]])

mean.n = colMeans(data.w[!y, 2:dim(data.w)[2]])
cov.n  =      cov(data.w[!y, 2:dim(data.w)[2]])

write.csv(sam, file.path('data', 'sim_parameters',
                         paste(d.set, 'sample_size.csv', sep = '_')),
          row.names = FALSE)
write.csv(mean.p, file.path('data', 'sim_parameters',
                            paste(d.set, 'Tclas', 'mean.csv', sep = '_')))
write.csv(cov.p,  file.path('data', 'sim_parameters',
                            paste(d.set, 'Tclas', 'cov.csv',  sep = '_')))
write.csv(mean.n, file.path('data', 'sim_parameters',
                            paste(d.set, 'Fclas', 'mean.csv', sep = '_')))
write.csv(cov.n,  file.path('data', 'sim_parameters',
                            paste(d.set, 'Fclas', 'cov.csv',  sep = '_')))
