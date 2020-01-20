
# Functions used in scripts. Mostly helper functions, plus lda_lw(), a 
# deterministic (as opposed to MASS:lda) implementation of LDA.
source('housekeeping_functions.R')

# Input from ./data/endometrial_cancer_data/
cat("\n---\n Generation of simulated data \n---\n")
source("sim_params.R")  # Calculate and store mean and covariance matrices for simulations
# Output to ./data/sim_parameters/ # Note: these files removed for repo for exceeding file-size restrictions but can be regenerated from inputs by running sim_params.R
source("sim_data.R")    # Simulate data from multivariate normal with above mean and covariance matrices.
# Output to ./data/sim_data

# Input from ./data/endometrial_cancer_data/
cat("\n---\n Analysis of endometrial cancer data \n---\n")
source("pca_lda.R")     # Perform PCA-LDA on endometrial cancer data (On various subsets of data, [Pi1] and [Pi2], etc.)
source("cca_lda.R")     # Perform CCA-LDA on endometrial cancer data (On various subsets of data, [Pi1] and [Pi2], etc.)
# Output to ./data/results/

# Input from ./data/sim_data
# Input from ./data/sim_parameters/
cat("\n---\n Analysis of simulated data \n---\n")
source("sim_pca_lda.R") # Perform PCA-LDA on simulated data (On various subsets of data, [Pi1] and [Pi2], etc.)
source("sim_cca_lda.R") # Perform CCA-LDA on simulated data (On various subsets of data, [Pi1] and [Pi2], etc.)
# Output to ./data/sim_results/

# Input from ./data/results
cat("\n---\n Summarise results on endometrial cancer data \n---\n")
source("summarise_results.R") # Condense and summarise results on endometrial cancer data
# Ouput to ./data/result_summaries

# Input from ./data/sim_results
cat("\n---\n Summarise results on simlated data \n---\n")
source("sim_summarise_results.R") # Condense and summarise results on simulated data
# Ouput to ./data/result_summaries

# Input from ./data/result_summaries
cat("\n---\n Produce figures for paper \n---\n")
source("paper_figures.R") # Produce figures
# Output to ../figures/


