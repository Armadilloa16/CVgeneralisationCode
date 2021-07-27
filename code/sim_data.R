
N = 200

sam = load_sim_param(TRUE, 'sample_size')

for (clas in c('p', 'n')) {
  cat(paste0("    class: ", clas, ", Simulation: "))
  
  # Load simulation parameters
  mean.tmp = load_sim_param(clas == 'p', 'mean')
  cov.tmp  = load_sim_param(clas == 'p', 'cov')
  
  for (sim in 1:N) {
    cat(paste(toString(sim), ' '))
    set.seed(as.numeric(paste0(toString(which(clas == c('p', 'n'))), 
                               toString(sim), '123')))
    
    # Simulate Data
    data.sim = mvrnorm(n = sam[, clas],
                       mu = mean.tmp,
                       Sigma = cov.tmp)
    
    # Save simulated data
    if (clas == 'p') {
      clas.str = '_Tclas'
    } else {
      clas.str = '_Fclas'
    }
    shift.str = '.csv'
    write.csv(data.sim, file.path("data", "sim_data",
                                  paste0('ETMA', clas.str,  
                                         '_', toString(sim), 
                                         shift.str)),
              row.names = FALSE)
  }
  cat('\n')
}

