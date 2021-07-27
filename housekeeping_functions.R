
library(geigen)



# Helper functions intended for managing application of scripts to different datasets

load_clin = function() {
  return(read.csv(file.path("data", "endometrial_cancer_data", "annotation", 
                            "clinical_variables","ETMA_clinical_variable.csv")))
}

load_patient_summary = function() {
    return(read.csv(file.path("data", "endometrial_cancer_data", 
                              "patient_summaries", 
                              "ETMA_logI_patient_averages.csv")))
}

reshape_data_matrix = function(pMeans, clin.df) {
  data.w = reshape(subset(pMeans, 
                          Patient %in% subset(clin.df, Suitable)$Patient),
                   timevar = "bin",
                   idvar = "Patient",
                   direction = "wide")
  
  # Replace missing values with zeroes.
  data.w[is.na(data.w)] = 0
  return(data.w)
}

get_sample_size = function(clin.df) {
  return(length(unique(subset(clin.df, Suitable)$Patient)))
}

get_clin_variable_name = function(d.set) {
  return("LNM")
}



# Helper functions intended to access stored data in systematically named files.

load_sim_param = function(clas, param, code_dir = '.') {
  d.set = "ETMA"
  shift.str = '.csv'
  if (clas) {
    clas.str = '_Tclas'
  } else {
    clas.str = '_Fclas'
  }
  param = paste0('_', param)
  if (param == '_mean') {
    return(read.csv(file.path(code_dir, "output", "simulation_parameters",
                              paste0(d.set, clas.str, param, shift.str)))[, 2])
  } else if (param == '_cov') {
    return(as.matrix(read.csv(file.path(code_dir,"output", "simulation_parameters",
                                        paste0(d.set, clas.str, param, shift.str)))[, -1]))
  } else if (param == '_sample_size') {
    return(read.csv(file.path(code_dir, "output", "simulation_parameters",
                                        paste0(d.set, param, '.csv'))))
  }
}

load_sim_data = function(clas, sim, code_dir = '.') {
  shift.str = '.csv'
  if (clas) {
    clas.str = '_Tclas'
  } else {
    clas.str = '_Fclas'
  }
  sim.str = paste0('_', toString(sim))
  return(read.csv(file.path(code_dir, "output", "simulated_data",
                            paste0('ETMA', clas.str, sim.str, shift.str))))
}



# Custom implementation of a deterministic (as opposed to MASS::lda) LDA 

train_lda = function(x, y, verbose = TRUE) {
  
  if (dim(x)[1] != length(y)) {
    stop('lda_lw: dimension mismatch')
  } else {
    d = dim(x)[2]
  }
  
  if (!is.logical(y)) {
    stop('lda_lw: requires y as logical vector')
  }
  
  # Calculate group (class) means, for group "TRUE" and group "FALSE"
  grp_means = data.frame(x[1:2, ], row.names = c("FALSE", "TRUE"))
  grp_means["FALSE", ] = colMeans(x[!y, ])
  grp_means["TRUE", ] = colMeans(x[y, ])
  
  # Between class covariance matrix
  B = cov(grp_means)
  
  # Sum of within class covariance matrices
  W = cov(x[!y, ]) + cov(x[y, ])
  
  # Eigendecomposition of W (only calculate eigenvalues)
  Weig = eigen(W, symmetric = TRUE, only.values = TRUE)
  
  if (any(abs(Weig$values) < 1e-14)) {
    
    # Rank of W
    r = min(which(Weig$values < 1e-14)[1]) - 1
    if (verbose) {
      cat(paste('\n lda_lw: r = ', r, 'when d =', d, 'next largest eigenvalue is', Weig$values[r], ''))
    }
    
    # Calculate full eigen-decomposition
    Weig = eigen(W, symmetric = TRUE)
    
    # Use the Moore-Penrose pseudo-inverse as a substitute for the inverse
    Winv = Weig$vectors[, 1:r] %*% diag(Weig$values[1:r]^-1) %*% t(Weig$vectors[, 1:r])
    
    # Multiply our pseudoinverse with B
    WinvB = Winv %*% B
    
    # Calculate the eigendecomposition of the product
    eig = eigen(WinvB)
    
    # Take the first eigenvector as our direction vector
    d = eig$vectors[, 1]
    
    # Small numeric errors can result in d being complex
    if (is.complex(d)) {
      if (verbose) {
        cat(paste('\n d complex with largest imaginary component', max(abs(Im(d))), ''))
      }
      d = Re(d)
    }
    
  } else {
    # Solve the generalised eigenvalue problem directly
    eig = geigen(B, W, symmetric = TRUE)
    d = eig$vectors[, which.max(eig$values)]
    
    # Normalise
    d = d / sqrt(sum(d^2))
  }
  
  # Project group means into direction d
  proj_grp_means = as.matrix(grp_means) %*% d
  
  # If neccessary multiply d by -1 so that the projected mean of group "TRUE" is 
  # greater than the projected mean of group "FALSE"
  if (proj_grp_means["FALSE", ] > proj_grp_means["TRUE", ]) {
    d = -1 * d
    proj_grp_means = as.matrix(grp_means) %*% d
  }
  
  # A new observation x will be classified into group "TRUE" if: 
  #   t(x) %*% d > cutoff, and into group "FALSE" if 
  #   t(x) %*% d <= cutoff.
  return(list(d = d, 
              cutoff = mean(proj_grp_means)))
  
}





