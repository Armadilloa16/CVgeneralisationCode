# CVgeneralisationCode
Code accompanying publication

Running script `workflow.R` should reproduce the analysis reported on in publication. This uses initial data in `./data/endometrial_cancer_data/`, and the remainder of the files in `./data/*` are intermediate files and results. Figures are output to `../figures/` by default so the final script run in `workflow.R`, `paper_figures.R`, will throw an error if this folder does not exist. All intermediate files and results other than the initial data in `./data/endometrial_cancer_data/` can be reproduced by running `workflow.R` but with the exception of `./data/sim_parameters/*` and `./data/result_summaries/SIM_Ecv6Pi1.csv` these intermediate files are also included in this repository for ease of access. The exceptions are ommitted only because of the default file-size limitation on GitHub.

Functions are included in `housekeeping_functions.R`, are loaded at the beggining of `workflow.R` and consist of helper functions and a deterministic implementation of LDA in `train_lda()`.

