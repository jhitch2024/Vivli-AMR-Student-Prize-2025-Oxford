# 0.Housekeeping ####

# 0.1 Clearing environment and graphics
rm(list = ls())
graphics.off()

# 0.2 Loading packages
library(tidyverse)
library(rstan)
library(bayesplot)
library(dplyr)
library(mvtnorm)
library(posterior)
library(tibble)
library(ggplot2)
library(forcats)
library(tableone)
library(cmdstanr)
library(AMR)
library(writexl)
library(arsenal)
library(furrr)
library(countrycode)
plan(multisession)

# 0.3 Setting working directory if desired
setwd("")

# 0.4 Reading in the final clean data for Southern Europe
load("vivli_data_final.RDS")
df_save <- vivli_data_final
rm(vivli_data_final)

# compiling the Bayesian models
mod_fixed <- cmdstan_model("mvprobit_fixed_final.stan") # non-hierarchical fixed effects model
mod_hier <- cmdstan_model("mvprobit_hierarchical_final.stan") # hierarchical model

# checking that the syntax is correct
mod_fixed$check_syntax()
mod_hier$check_syntax()

# new data table with only model variables
df <- df_save %>%
  rename(icu = icu_bin)
df <- df %>%
  dplyr::select(country, region, year, gender, age, icu, 
         Meropenem_bin, 
         Piperacillin.tazobactam_bin,
         Ampicillin_bin, 
         Tigecycline_bin,
         Gentamicin_bin, 
         Amikacin_bin, 
         Ceftazidime_bin,
         Levofloxacin_bin, 
         Colistin_bin, 
         Amoxycillin.clavulanate_bin)
summary(df)

# converting metadata variables to factors
df <- df %>%
  mutate(across(c(country, region, year, gender, age, icu), as.factor))

# splitting the data into training (2018-2022) and test sets (2023)
df_south_train <- df %>%
  filter(year %in% c(2018:2022))
df_south_test <- df %>%
  filter(year == 2023)

# creating a duplicate of the test set for fitting the fixed effects model
df_south_test_fixed <- df_south_test

# dropping unused levels from country to ensure they don't appear in model matrix
df_south_train$country <- droplevels(df_south_train$country)
df_south_test$country <- droplevels(df_south_test$country)

# specifying training design matrices
# The design matrix for the fixed effect model will differ from the hierarchical 
# model in that it will include dummies for country whereas the design matrix
X_south_fixed <- model.matrix(~ age + gender + icu + country, data = df_south_train)[,-1]
X_south_hier <- model.matrix(~ age + gender + icu, data = df_south_train)[,-1]

# specifying test design matrices
X_south_fixed_test <- model.matrix(~ age + gender + icu + country, data = df_south_test)[,-1]
X_south_hier_test <- model.matrix(~ age + gender + icu, data = df_south_test)[,-1]

# specifying outcome matrix (same for fixed and hierarhical models)
Y_south <- as.matrix(df_south_train[, c("Meropenem_bin", 
                                        "Piperacillin.tazobactam_bin",
                                        "Ampicillin_bin", 
                                        "Tigecycline_bin",
                                        "Gentamicin_bin", 
                                        "Amikacin_bin", 
                                        "Ceftazidime_bin",
                                        "Levofloxacin_bin", 
                                        "Colistin_bin", 
                                        "Amoxycillin.clavulanate_bin")])


# declaring the data to fit the STAN fixed effect model
stan_data_south_fixed <- list(
  N = nrow(df_south_train),
  P = ncol(X_south_fixed),       
  D = ncol(Y_south),             
  X = X_south_fixed,          
  Y = Y_south,        
  N_test =  nrow(df_south_test),
  X_test = X_south_fixed_test
)

# fitting the STAN fixed effects model
fit_fixed_south <- mod_fixed$sample(
  data = stan_data_south_fixed,
  chains = 4, # change this to 4
  parallel_chains = 4, # change this to 4
  iter_warmup = 3000, # change this to 2500
  iter_sampling = 3000, # change this to 2500
  seed = 1954,
  adapt_delta = 0.95)

# running diagnostics
fit_fixed_south$cmdstan_diagnose()

# saving the model
fit_fixed_south$save_object(file = "model_outputs/fit_fixed_south.RDS")

params_problematic <- c(
  "L[9,2]", "L[4,3]", "L[9,3]", "L[10,3]",
  "L[4,4]", "L[9,4]"
)

draws_subset <- fit_fixed_south$draws(variables = params_problematic, format = "df")
draws_array <- fit_fixed_south$draws(variables = params_problematic, format = "draws_array")

alpha_array <- fit_hier_south$draws("sigma_alpha")
summary(alpha_array)

mcmc_trace(draws_array)

fit_fixed_south$draws("L[9,2]", format = "draws_array")
fit_fixed_south$draws("L[9,3]", format = "draws_array")

summary_df <- fit_fixed_south$summary(variables = params_problematic)
summary_df[, c("variable", "rhat")]

fit_fixed_south$summary(variables = params_problematic, 
                        summary = c("mean", "sd", "ess_bulk", "ess_tail", "rhat"))

# declaring the data to fit the STAN hierarchical model
# firstly we have to convert country into an integer factor variable for the
# hierarchical model to work
df_south_train <- df_south_train %>%
  mutate(country_id = as.integer(factor(country)))
df_south_test <- df_south_test %>%
  mutate(country_id_test = as.integer(factor(country, levels = levels(factor(df_south_train$country))))) # this ensures they match in the training and test sets

# creating a duplicate of the test set for fitting the hierarchical model
df_south_test_hier <- df_south_test

# declaring the data to fit the STAN hierarchical model
stan_data_south_hier <- list(
  N = nrow(df_south_train),
  P = ncol(X_south_hier),     
  D = ncol(Y_south),             
  X = X_south_hier,          
  Y = Y_south,         
  C = length(unique(df_south_train$country)),
  country_id = as.integer(df_south_train$country),
  N_test =  nrow(df_south_test),
  X_test = X_south_hier_test,
  country_id_test = as.integer(df_south_test$country) # Ensure 1-based indexing
)

# fitting the STAN hierarchical model
fit_hier_south <- mod_hier$sample(
  data = stan_data_south_hier,
  chains = 4, # change this to 4
  parallel_chains = 4, # change this to 4
  iter_warmup = 3000, # change this to 3000
  iter_sampling = 3000, # change this to 3000
  seed = 1954,
  adapt_delta = 0.95)

# running diagnostics
fit_hier_south$cmdstan_diagnose()

params_problematic <- c(
  "L[4,2]", "L[9,2]", "L[4,3]",
  "L[9,3]", "L[10,3]"
)

draws_subset_hier <- fit_hier_south$draws(variables = params_problematic, format = "df")
draws_array_hier <- fit_hier_south$draws(variables = params_problematic, format = "draws_array")

mcmc_trace(draws_array_hier)


summary_df <- fit_hier_south$summary(variables = params_problematic)
summary_df[, c("variable", "rhat")]

fit_hier_south$summary(variables = params_problematic, 
                        summary = c("mean", "sd", "ess_bulk", "ess_tail", "rhat"))


if (!dir.exists("model_outputs")) {
  dir.create("model_outputs")
}


######################### DRAWS #############

draws_pred_prob_df_fixed <- fit_fixed_south$draws("pred_prob_test", format = "df")
y_test_rep_fixed <- fit_fixed_south$draws("Y_test_rep", format = "df")
z_test_fixed <- fit_fixed_south$draws("Z_test", format = "df")

alpha_fixed <- fit_fixed_south$draws("alpha", format = "df")
beta_fixed <-  fit_fixed_south$draws("beta", format = "df")
omega_fixed <-  fit_fixed_south$draws("Omega", format = "df")

saveRDS(alpha_fixed, file = "model_outputs/alpha_fixed.RDS")
saveRDS(beta_fixed, file = "model_outputs/beta_fixed.RDS")
saveRDS(omega_fixed, file = "model_outputs/omega_fixed.RDS")

alpha_hier <- fit_hier_south$draws("alpha", format = "df")
beta_hier <-  fit_hier_south$draws("beta", format = "df")
omega_hier <-  fit_hier_south$draws("Omega", format = "df")

saveRDS(alpha_hier, file = "model_outputs/alpha_hier.RDS")
saveRDS(beta_hier, file = "model_outputs/beta_hier.RDS")
saveRDS(omega_hier, file = "model_outputs/omega_hier.RDS")


saveRDS(y_test_rep_fixed, file = "model_outputs/y_test_rep_fixed_v3.RDS")
saveRDS(draws_pred_prob_df_fixed, file = "model_outputs/draws_pred_prob_df_fixed_v3.RDS")
saveRDS(z_test_fixed, file = "model_outputs/z_test_fixed.RDS")

saveRDS(alpha_fixed, file = "model_outputs/alpha_fixedRDS")



draws_pred_prob_df_hier <- fit_hier_south$draws("pred_prob_test", format = "df")
y_test_rep_hier <- fit_hier_south$draws("Y_test_rep", format = "df")

z_test_hier <- fit_hier_south$draws("Z_test", format = "df")



saveRDS(y_test_rep_hier, file = "model_outputs/y_test_rep_hier_v3.RDS")
saveRDS(draws_pred_prob_df_hier, file = "model_outputs/draws_pred_prob_df_hier_v3.RDS")



saveRDS(joint_data_hier_model, file = "model_outputs/joint_data_hier_model.RDS")
saveRDS(joint_data_fixed_model, file = "model_outputs/joint_data_fixed_model.RDS")
saveRDS(z_test_hier, file = "model_outputs/z_test_hier.RDS")


draws_pred_prob_df_hier <- readRDS("model_outputs/draws_pred_prob_df_hier_v3.RDS")


##################### Probabilities estimation ##########################

# Create a list with the antibiotics name

ab_names <- c("Meropenem_bin",
              "Piperacillin.tazobactam_bin",
              "Ampicillin_bin",
              "Tigecycline_bin",
              "Gentamicin_bin",
              "Amikacin_bin",
              "Ceftazidime_bin",
              "Levofloxacin_bin",
              "Colistin_bin",
              "Amoxycillin.clavulanate_bin")


prepare_joint_dataset <- function(fit, df, ab_names) {
  df <- df %>% mutate(patient_id = row_number())
  N <- nrow(df)
  D <- length(ab_names)
  
  # Extract Z_test draws (latent variables)
  zrep_raw <- fit$draws("Z_test", format = "df")
  
  # Pivot to long format
  zrep_df <- zrep_raw %>%
    pivot_longer(cols = starts_with("Z_test["), 
                 names_to = "var", values_to = "value") %>% 
    mutate(
      var = str_remove_all(var, "Z_test\\[|\\]"),
      patient_id = as.integer(str_extract(var, "^[0-9]+")),
      ab_index = as.integer(str_extract(var, "(?<=,)[0-9]+")),
      ab = ab_names[ab_index]
    ) %>%
    dplyr::select(.draw, patient_id, ab, value)
  
  # Pivot wide: one row per patient-draw with antibiotics as columns
  zrep_wide <- zrep_df %>%
    pivot_wider(names_from = ab, values_from = value)
  
  # Add covariates
  df_covars <- df %>%
    dplyr::select(patient_id, age, gender, country, icu, year)
  
  final_df <- zrep_wide %>%
    left_join(df_covars, by = "patient_id")
  
  return(final_df)
}

joint_data_hier_model <- prepare_joint_dataset(
  fit_hier_south, 
  df_south_test_hier, 
  ab_names
)

estimate_conditional_probs_patient_first_resistant <- function(joint_data, ab_names) {
  # Generate all possible antibiotic pairs (A ≠ B)
  combos <- expand.grid(ab_A = ab_names, ab_B = ab_names, stringsAsFactors = FALSE) %>%
    dplyr::filter(ab_A != ab_B)
  
  results <- purrr::map_dfr(1:nrow(combos), function(i) {
    ab_A <- combos$ab_A[i]
    ab_B <- combos$ab_B[i]
    
    joint_data %>%
      dplyr::mutate(
        # Convert latent Z into probability of resistance (Probit link)
        prob_A = pnorm(!!rlang::sym(ab_A)),  
        prob_B = pnorm(!!rlang::sym(ab_B))   
      ) %>%
      # Condition on likely resistance to A
      dplyr::filter(prob_A >= 0.5) %>%  
      dplyr::group_by(.draw, age, gender, country, icu, year) %>%
      dplyr::summarise(
        cond_prop = mean(prob_B),  # Expected P(resistant to B)
        .groups = "drop"
      ) %>%
      dplyr::group_by(age, gender, country, icu, year) %>%
      dplyr::summarise(
        cond_ab_given = paste0("P(", ab_B, "=1 | ", ab_A, "=1)"),
        ab1 = ab_A,
        ab2 = ab_B,
        prob_cond = mean(cond_prop),
        ci_lower = quantile(cond_prop, 0.025),
        ci_upper = quantile(cond_prop, 0.975),
        .groups = "drop"
      )
  })
  
  return(results)
}

cond_results_hier_model <- estimate_conditional_probs_patient_first_resistant(joint_data_hier_model, ab_names)

saveRDS(cond_results_hier_model, file = "model_outputs/cond_results_hier_model.RDS")


estimate_conditional_probs_patient_first_resistant_to_susceptible <- function(joint_data, ab_names) {
  # Generate all possible antibiotic pairs (A ≠ B)
  combos <- expand.grid(ab_A = ab_names, ab_B = ab_names, stringsAsFactors = FALSE) %>%
    dplyr::filter(ab_A != ab_B)
  
  results <- purrr::map_dfr(1:nrow(combos), function(i) {
    ab_A <- combos$ab_A[i]
    ab_B <- combos$ab_B[i]
    
    joint_data %>%
      dplyr::mutate(
        prob_A = pnorm(!!rlang::sym(ab_A)),  # P(resistant to A)
        prob_B = pnorm(!!rlang::sym(ab_B))   # P(resistant to B)
      ) %>%
      dplyr::filter(prob_A >= 0.5) %>%  
      dplyr::group_by(.draw, age, gender, country, icu, year) %>%
      dplyr::summarise(
        cond_prop = 1 - mean(prob_B),  # Expected P(resistant to B)
        .groups = "drop"
      ) %>%
      dplyr::group_by(age, gender, country, icu, year) %>%
      dplyr::summarise(
        cond_ab_given = paste0("P(", ab_B, "=0 | ", ab_A, "=1)"),
        ab1 = ab_A,
        ab2 = ab_B,
        prob_cond = mean(cond_prop),
        ci_lower = quantile(cond_prop, 0.025),
        ci_upper = quantile(cond_prop, 0.975),
        .groups = "drop"
      )
  })
  
  return(results)
}


cond_results_hier_model_s <- estimate_conditional_probs_patient_first_resistant_to_susceptible(joint_data_hier_model, ab_names)
  
saveRDS(cond_results_hier_model_s, file = "model_outputs/cond_results_hier_model_susceptible.RDS")

# Conditional probabilities by country 

estimate_conditional_probs_by_country <- function(joint_data, ab_names) {
  # Generate all possible antibiotic pairs (A ≠ B)
  combos <- expand.grid(ab_A = ab_names, ab_B = ab_names, stringsAsFactors = FALSE) %>%
    dplyr::filter(ab_A != ab_B)
  
  # Iterate over each pair to estimate P(B = 0 | A = 1) at country level
  results <- purrr::map_dfr(1:nrow(combos), function(i) {
    ab_A <- combos$ab_A[i]
    ab_B <- combos$ab_B[i]
    
    joint_data %>%
      dplyr::mutate(prob_A = pnorm(!!rlang::sym(ab_A)),
                    prob_B = pnorm(!!rlang::sym(ab_B))) %>%
      dplyr::filter(prob_A >= 0.5) %>%  # Condition on A being resistant
      dplyr::group_by(.draw, country) %>%
      dplyr::summarise(
        cond_prop = 1 - mean(prob_B),   # P(B = 0 | A = 1)
        .groups = "drop"
      ) %>%
      dplyr::group_by(country) %>%
      dplyr::summarise(
        cond_ab_given = paste0("P(", ab_B, "=0 | ", ab_A, "=1)"),
        ab1 = ab_A,
        ab2 = ab_B,
        prob_cond = mean(cond_prop),
        ci_lower = quantile(cond_prop, 0.025),
        ci_upper = quantile(cond_prop, 0.975),
        .groups = "drop"
      )
  })
  
  return(results)
}


cond_results_hier_model_by_country <- estimate_conditional_probs_by_country(joint_data_hier_model, ab_names)


cond_results_clean <- cond_results_hier_model %>%
  mutate(
    profile_id = paste(age, gender, country, icu, year, sep = "_"),
    var_base = paste0(ab2, "_given_", ab1)
  )

wide_mean <- cond_results_clean %>%
  dplyr::select(profile_id, var_base, prob_cond) %>%
  pivot_wider(names_from = var_base, values_from = prob_cond, names_prefix = "mean_")

wide_lower <- cond_results_clean %>%
  dplyr::select(profile_id, var_base, ci_lower) %>%
  pivot_wider(names_from = var_base, values_from = ci_lower, names_prefix = "lower_")

wide_upper <- cond_results_clean %>%
  dplyr::select(profile_id, var_base, ci_upper) %>%
  pivot_wider(names_from = var_base, values_from = ci_upper, names_prefix = "upper_")

profile_info <- cond_results_clean %>%
  dplyr::select(profile_id, age, gender, country, icu, year) %>%
  distinct()

final_conditional_results <- profile_info %>%
  left_join(wide_mean,  by = "profile_id") %>%
  left_join(wide_lower, by = "profile_id") %>%
  left_join(wide_upper, by = "profile_id") %>%
  dplyr::select(-profile_id)

saveRDS(final_conditional_results, file = "model_outputs/cond_results_hier_model_resistant_susceptible.RDS")

res_patient_level_hier <- estimate_conditional_probs_per_patient(joint_data_hier_model, ab_names)

res_summary_hier <- res_patient_level_hier %>%
  group_by(patient_id, ab1, ab2) %>%
  summarise(
    prob_cond = mean(prob_cond),
    ci_lower = quantile(prob_cond, 0.025),
    ci_upper = quantile(prob_cond, 0.975),
    .groups = "drop"
  )

res_summary_hier <- res_summary_hier %>%
  mutate(
    col_prob  = paste0("P(", ab2, "=1 | ", ab1, "=1)"),
    col_lower = paste0("P(", ab2, "=1 | ", ab1, "=1)_lower"),
    col_upper = paste0("P(", ab2, "=1 | ", ab1, "=1)_upper")
  )

res_wide_prob_hier <- res_summary_hier %>%
  select(patient_id, col_prob, prob_cond) %>%
  pivot_wider(names_from = col_prob, values_from = prob_cond)

res_wide_lower_hier <- res_summary_hier %>%
  select(patient_id, col_lower, ci_lower) %>%
  pivot_wider(names_from = col_lower, values_from = ci_lower)

res_wide_upper_hier <- res_summary_hier %>%
  select(patient_id, col_upper, ci_upper) %>%
  pivot_wider(names_from = col_upper, values_from = ci_upper)

res_joint_hier <- reduce(list(res_wide_prob_hier, res_wide_lower_hier, res_wide_upper_hier), full_join, by = "patient_id")

df_final_hier <- df_south_test_hier %>%
  mutate(patient_id = row_number()) %>%
  left_join(res_joint_hier, by = "patient_id")



##### Visualization #######


cond_results_hier_model_plot <- cond_results_hier_model %>%
  mutate(
    patient_profile = paste(age, gender, ifelse(icu == 1, "ICU", "No ICU"), sep = " | ")
  )

cond_results_hier_model_plot <- cond_results_hier_model_plot %>%
  dplyr::mutate(facet_label = paste(country, patient_profile, sep = " - "))

cond_results_hier_model_plot <- cond_results_hier_model_plot %>% 
  mutate(
    ab1_clean = str_remove(ab1, "_bin$"),
    ab2_clean = str_remove(ab2, "_bin$")
  ) %>%   mutate(
    ab1_abbrev = recode(ab1_clean,
                        "Piperacillin.tazobactam" = "TZP",
                        "Ampicillin" = "AMP",
                        "Tigecycline" = "TGC",
                        "Gentamicin" = "GEN",
                        "Amikacin" = "AMK",
                        "Ceftazidime" = "CAZ",
                        "Levofloxacin" = "LVX",
                        "Colistin" = "CST",
                        "Amoxycillin.clavulanate" = "AMC",
                        "Meropenem" = "MEM"
    ),
    ab2_abbrev = recode(ab2_clean,
                        "Piperacillin.tazobactam" = "TZP",
                        "Ampicillin" = "AMP",
                        "Tigecycline" = "TGC",
                        "Gentamicin" = "GEN",
                        "Amikacin" = "AMK",
                        "Ceftazidime" = "CAZ",
                        "Levofloxacin" = "LVX",
                        "Colistin" = "CST",
                        "Amoxycillin.clavulanate" = "AMC",
                        "Meropenem" = "MEM"
    )
  )



p1 <- ggplot(cond_results_hier_model_plot, aes(x = ab1_abbrev, y = ab2_abbrev, fill = prob_cond)) +
  geom_tile(color = "white") +
  facet_wrap(~ facet_label) +
  scale_fill_viridis_c(name = "P(B = 0 | A = 1)", limits = c(0,1)) +
  theme_minimal(base_size = 10) +
  labs(
    title = "Estimated Conditional Probabilities of Susceptibility - Hierarchical model",
    subtitle = "Derived from Multivariate Probit Model",
    x = "Antibiotic A (conditioned on A = 1)",
    y = "Antibiotic B (probability of B = 0)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave(
  "model_outputs/conditional_probs_by_profile_hierarchial_model_V2.png",
  plot = p1,
  width = 24,   
  height = 16,  
  dpi = 300,
  units = "in"
)


cond_results_hier_model_by_country <-  cond_results_hier_model_by_country %>% 
  mutate(prob_cond = prob_cond * 100,
         ci_lower = ci_lower * 100,
         ci_upper = ci_upper * 100)

cond_results_hier_model_by_country <- cond_results_hier_model_by_country %>% 
  mutate(
    ab1_clean = str_remove(ab1, "_bin$"),
    ab2_clean = str_remove(ab2, "_bin$")
  ) %>%   mutate(
    ab1_abbrev = recode(ab1_clean,
                        "Piperacillin.tazobactam" = "TZP",
                        "Ampicillin" = "AMP",
                        "Tigecycline" = "TGC",
                        "Gentamicin" = "GEN",
                        "Amikacin" = "AMK",
                        "Ceftazidime" = "CAZ",
                        "Levofloxacin" = "LVX",
                        "Colistin" = "CST",
                        "Amoxycillin.clavulanate" = "AMC",
                        "Meropenem" = "MEM"
    ),
    ab2_abbrev = recode(ab2_clean,
                        "Piperacillin.tazobactam" = "TZP",
                        "Ampicillin" = "AMP",
                        "Tigecycline" = "TGC",
                        "Gentamicin" = "GEN",
                        "Amikacin" = "AMK",
                        "Ceftazidime" = "CAZ",
                        "Levofloxacin" = "LVX",
                        "Colistin" = "CST",
                        "Amoxycillin.clavulanate" = "AMC",
                        "Meropenem" = "MEM"
    )
  )




p2 <- ggplot(cond_results_hier_model_by_country, aes(x =ab1_abbrev, y = ab2_abbrev, fill = prob_cond)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n[%.2f, %.2f]", prob_cond, ci_lower, ci_upper)), check_overlap = F, colour = 'white',size = 2.5) +
  facet_wrap(~ country, ncol = 3) +
  scale_fill_viridis_c(name = "P(B = 0 | A = 1)", limits = c(0, 100)) +
  labs(
    title = "Conditional Susceptibility Probabilities by Country - Multivariate  Hiearchical Probit Model",
    subtitle = "Mean and 95% CrI",
    x = "Resistant to Antibiotic A (A = 1)",
    y = "Susceptible to Antibiotic B (P(B = 0))"
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 14)
  )

ggsave(
  "model_outputs/conditional_probs_by_country_hierarchial_model.png",
  plot = p2,
  width = 24,   
  height = 16,  
  dpi = 300,
  units = "in"
)


################ Processing Z draws for presentation #############


prepare_z_dataset <- function(z_draws, df, ab_names, n_draws = NULL) {
  # Ensure unique ID per patient
  df <- df %>% mutate(patient_id = row_number())
  N <- nrow(df)
  D <- length(ab_names)
  
  # Filter draws if requested
  if (!is.null(n_draws)) {
    z_draws <- z_draws %>% slice(1:n_draws)
  }
  
  # Pivot to long format
  z_long <- z_draws %>%
    pivot_longer(cols = starts_with("Z_test["), names_to = "var", values_to = "z_value") %>% 
    mutate(
      var = str_remove_all(var, "Z_test\\[|\\]"),
      patient_id = as.integer(str_extract(var, "^[0-9]+")),
      ab_index   = as.integer(str_extract(var, "(?<=,)[0-9]+")),
      antibiotic = ab_names[ab_index]
    ) %>%
    dplyr::select(.draw, patient_id, antibiotic, z_value)
  
  # Add covariates
  df_covars <- df %>%
    dplyr::select(patient_id, age, gender, country, icu, year)
  
  final_df <- z_long %>%
    left_join(df_covars, by = "patient_id")
  
  return(final_df)
}

df_z_long <- prepare_z_dataset(z_test_hier, df_south_test, ab_names, n_draws = 12000)
cond_results_hier_model_z <- estimate_conditional_probs_patient_first(df_z_long, ab_names)

estimate_conditional_probs_by_profile <- function(df_z_long, ab_A, ab_B, countries = NULL, n_draws = NULL) {
  
  df_filtered <- df_z_long %>%
    filter(antibiotic %in% c(ab_A, ab_B)) %>%
    { if (!is.null(countries)) filter(., country %in% countries) else . } %>%
    { if (!is.null(n_draws)) filter(., .draw %in% sample(unique(.draw), n_draws)) else . }
  
  df_wide <- df_filtered %>%
    tidyr::pivot_wider(
      id_cols = c(.draw, patient_id, age, gender, icu, year, country),
      names_from = antibiotic,
      values_from = z_value
    )
  
  df_result <- df_wide %>%
    dplyr::mutate(
      prob_A = pnorm(.data[[ab_A]]),
      prob_B = pnorm(.data[[ab_B]]),
      cond_prob = ifelse(prob_A >= 0.5, 1 - prob_B, NA_real_), # Probability of B susceptible given A resistant
      cond_ab = paste0("P(", ab_B, "=0 | ", ab_A, "=1)")
    ) %>%
    dplyr::filter(!is.na(cond_prob))  # Keep only those that meet the condition
  
  return(df_result)
}

# Estimate conditional probabilities by profile
df_by_profile <- estimate_conditional_probs_by_profile(
  df_z_long,
  ab_A = "Piperacillin.tazobactam_bin",
  ab_B = "Meropenem_bin",
  countries = c("Italy", "Spain", "Croatia",
                "Greece", "Portugal", "Slovenia"),
  n_draws = 12000
)

# Recode ICU indicator
df_by_profile <- df_by_profile %>%
  mutate(
    icu = recode(icu,
                 "0" = "No ICU",
                 "1" = "ICU"
    )
  )

# Summarize posterior distributions by patient profile
df_summary <- df_by_profile %>%
  group_by(country, age, gender, icu) %>%
  summarise(
    median  = median(cond_prob),
    lower = quantile(cond_prob, 0.025),
    upper = quantile(cond_prob, 0.975),
    .groups = "drop"
  )

# Create profile label for y-axis
df_summary <- df_summary %>%
  mutate(profile = paste(country, age, gender, icu, sep = " | "))

# Plot: point = posterior median, line = 95% CrI
ggplot(df_summary, aes(x = median, y = reorder(profile, median),
                       colour = country)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Probability (median, 95% CrI)",
    y = "Country | Age | Gender | ICU",
    title = "Posterior median and 95% credible intervals by patient profile",
    colour = "Country"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

# Maximum line height for vertical annotations in density plots
line_height <- 2.5

df_by_profile %>%
  filter(country == "Italy") %>%
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  
  # Solid line = posterior median
  geom_segment(
    data = df_summary %>% filter(country == "Italy"),
    aes(x = median, xend = median, y = 0, yend = line_height),
    color = "black", size = 1
  ) +
  
  # Dashed lines = 95% CrI
  geom_segment(
    data = df_summary %>% filter(country == "Italy"),
    aes(x = lower, xend = lower, y = 0, yend = line_height),
    color = "black", linetype = "dashed", size = 0.8
  ) +
  geom_segment(
    data = df_summary %>% filter(country == "Italy"),
    aes(x = upper, xend = upper, y = 0, yend = line_height),
    color = "black", linetype = "dashed", size = 0.8
  ) +
  
  # Text label with median and CrI
  geom_text(
    data = df_summary %>% filter(country == "Italy"),
    aes(
      x = median,
      y = line_height + 0.2,
      label = sprintf("%.2f (%.2f–%.2f)", median, lower, upper)
    ),
    size = 3.2,
    hjust = 0.5
  ) +
  facet_grid(rows = vars(age, gender), cols = vars(icu)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacillin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()


df_by_profile %>%
  filter(country == "Spain") %>% 
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_grid(rows = vars(age, gender), cols = vars(icu)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacilin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()

df_by_profile %>%
  filter(country == "Croatia") %>% 
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_grid(rows = vars(age, gender), cols = vars(icu)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacilin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()

df_by_profile %>%
  filter(country == "Greece") %>% 
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_grid(rows = vars(age, gender), cols = vars(icu)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacilin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()

df_by_profile %>%
  filter(country == "Portugal") %>% 
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_grid(rows = vars(age, gender), cols = vars(icu)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacilin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()

df_by_profile %>%
  filter(country == "Slovenia") %>% 
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  facet_grid(rows = vars(age, gender), cols = vars(icu)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacilin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()



line_height <- 11.4  

df_by_profile %>%
  filter(country == "Slovenia") %>%
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
    geom_segment(
    data = df_summary %>% filter(country == "Slovenia"),
    aes(x = median, xend = median, y = 0, yend = line_height),
    color = "black", size = 1
  ) +
    geom_segment(
    data = df_summary %>% filter(country == "Slovenia"),
    aes(x = lower, xend = lower, y = 0, yend = line_height),
    color = "black", linetype = "dashed", size = 0.8
  ) +
  geom_segment(
    data = df_summary %>% filter(country == "Slovenia"),
    aes(x = upper, xend = upper, y = 0, yend = line_height),
    color = "black", linetype = "dashed", size = 0.8
  ) +
    geom_text(
    data = df_summary %>% filter(country == "Slovenia"),
    aes(
      x = median,
      y = line_height + 0.2,
      label = sprintf("%.2f (%.2f–%.2f)", median, lower, upper)
    ),
    size = 3.2,
    hjust = 0.5
  ) +
  
  facet_grid(rows = vars(age, gender), cols = vars(icu)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacillin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()


df_by_profile %>%
  filter(country == "Slovenia" & gender == "Female" & age == "19 to 64 Years") %>%
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
    geom_segment(
    data = df_summary %>% filter(country == "Slovenia" & gender == "Female" & age == "19 to 64 Years"),
    aes(x = median, xend = median, y = 0, yend = line_height),
    color = "black", size = 1
  ) +
    geom_text(
    data = df_summary %>% filter(country == "Slovenia" & gender == "Female" & age == "19 to 64 Years"),
    aes(
      x = median,
      y = line_height + 0.5,
      label = sprintf("%.2f (%.2f–%.2f)", median, lower, upper)
    ),
    size = 3.2,
    hjust = 0.5
  ) +
  
  facet_grid(cols = vars(icu)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacillin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()



df_by_profile %>%
  filter(country == "Spain" & age == "65 and Over" & icu == "No ICU") %>%
  ggplot(aes(x = cond_prob)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
    geom_segment(
    data = df_summary %>% filter(country == "Spain" & age == "65 and Over" & icu == "No ICU"),
    aes(x = median, xend = median, y = 0, yend = line_height),
    color = "black", size = 1
  ) +
    geom_text(
    data = df_summary %>% filter(country == "Spain" & age == "65 and Over" & icu == "No ICU"),
    aes(
      x = median,
      y = line_height + 0.5,
      label = sprintf("%.2f (%.2f–%.2f)", median, lower, upper)
    ),
    size = 3.2,
    hjust = 0.5
  ) +
  
  facet_grid(rows = vars(gender)) +
  labs(
    title = "Conditional probability: P(Meropenem = 0 | Piperacillin/tazobactam = 1)",
    x = "Conditional probability",
    y = "Density"
  ) +
  theme_minimal()
