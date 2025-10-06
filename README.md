# 2025 Vivli AMR Surveillance Data Challenge
This repository contains the code to reproduce our submission for the 2025 Vivli AMR Surveillance Data Challenge.

R_script_data_processing.R cleans the raw ATLAS data and derives our final analytical sample of gram-negative blood isolates collected across Southern Europe 2018-2023.

mvprobit_fixed_final.stan is the STAN code behind the fixed effects (FE) Bayesian multivariate probit model

mvprobit_fixed_hierarchical.stan is the STAN code behind the hierarchical Bayesian multivariate probit model (HB-MVP)

cond_probs_mvnorm.R takes the posterior draws from the two Bayesian regression models and derives predictions for susceptibility conditional on observed resistance/susceptibility to each index antibiotic, for all isolates in our analytical sample (these predictions are used to compare model performance)

