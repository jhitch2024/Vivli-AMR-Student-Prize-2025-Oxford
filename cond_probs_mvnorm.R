# setting working directory
rm(list = ls())
graphics.off()

# setting working directory if required
setwd("")

# reading in design matrices for the two models (no intercept)
X_south_fixed_test <- readRDS("X_south_fixed_test.RDS")
X_south_hier_test <- readRDS("X_south_hier_test.RDS")

# reading in "ground truth" outcome values
Y_south_test <- readRDS("Y_south_test.RDS")
Y_south_test <- Y_south_test[, colnames(Y_south_test) != "id"]

# combining covariates and outcomes separately for each model
X_Y_fixed <- cbind(X_south_fixed_test, Y_south_test)
X_Y_hier <- cbind(X_south_hier_test, Y_south_test)

# now have two dataframes with outcome and covariate values for 2023, one set 
# for each model

# focusing on the fixed effects model to start
# reading in parameters from fixed model
alpha_fixed <- readRDS("model_outputs/alpha_fixed.RDS")
beta_fixed <- readRDS("model_outputs/beta_fixed.RDS")
omega_fixed <- readRDS("model_outputs/omega_fixed.RDS")

alpha_hier <- readRDS("model_outputs/alpha_hier.RDS")
beta_hier <- readRDS("model_outputs/beta_hier.RDS")
omega_hier <- readRDS("model_outputs/omega_hier.RDS")

# dropping draw, iteration and chain from each because won't need this
alpha_fixed <- alpha_fixed %>%
  select(!c(".draw", ".iteration", ".chain"))
alpha_fixed_means <- colMeans(alpha_fixed)

beta_fixed <- beta_fixed %>%
  select(!c(".draw", ".iteration", ".chain"))
beta_fixed_means <- colMeans(beta_fixed)

omega_fixed <- omega_fixed %>%
  select(!c(".draw", ".iteration", ".chain"))
omega_fixed_means <- colMeans(omega_fixed)
params_fixed <- cbind(alpha_fixed, beta_fixed, omega_fixed)


alpha_hier <- alpha_hier %>%
  select(!c(".draw", ".iteration", ".chain"))
alpha_hier_means <- colMeans(alpha_hier)

beta_hier <- beta_hier %>%
  select(!c(".draw", ".iteration", ".chain"))
beta_hier_means <- colMeans(beta_hier)

omega_hier <- omega_hier %>%
  select(!c(".draw", ".iteration", ".chain"))
omega_hier_means <- colMeans(omega_hier)
params_hier <- cbind(alpha_hier, beta_hier, omega_hier)

params_fixed <- list(alpha_fixed, beta_fixed, omega_fixed)
params_fixed_list <- lapply(params_fixed, as.matrix)

params_hier <- list(alpha_hier, beta_hier, omega_hier)
params_hier_list <- lapply(params_hier, as.matrix)

# density plots for intercepts
mcmc_areas(alpha_fixed) # one per antibiotic
mcmc_areas(alpha_hier) # one per antibiotic per country

# density plots for beta coefficients
# meropenem as an example
head(beta_hier)
beta_fixed_mem <- beta_fixed %>% select(matches("beta\\[\\d+,1\\]"))
beta_hier_mem <- beta_hier %>% select(matches("beta\\[\\d+,1\\]"))


pars_fixed <- c("beta[1,1]", "beta[2,1]", "beta[3,1]", "beta[4,1]", "beta[5,1]", "beta[6,1]", "beta[7,1]")
labels_fixed <- c("Age (>65)", "Male", "ICU", "Greece vs Croatia", "Italy vs Croatia", "Portugal vs Croatia", "Slovenia vs Croatia", "Spain vs Croatia")

# Posterior density plots
mcmc_areas(beta_fixed_mem, 
           pars = pars_fixed, 
           prob = 0.8) +
  scale_y_discrete(labels = labels_fixed) +
  ggtitle("Posterior densities of beta coefficients for meropenem: fixed effects model")


pars_hier <- c("beta[1,1]", "beta[2,1]", "beta[3,1]")
labels_hier <- c("Age (>65)", "Male", "ICU")

mcmc_areas(beta_hier_mem, 
           pars = pars_hier, 
           prob = 0.8) +
  scale_y_discrete(labels = labels_hier) +
  ggtitle("Posterior densities of beta coefficients for meropenem: hierarchical model")



# correlation matrices
omega_fixed_means <- colMeans(omega_fixed)
omega_hier_means <- colMeans(omega_hier)

library(stringr)

# Extract i,j indices from names like Omega[2,1]
index_pairs <- str_match(names(omega_hier_means), "Omega\\[(\\d+),(\\d+)\\]")
i_vals <- as.integer(index_pairs[, 2])
j_vals <- as.integer(index_pairs[, 3])

# Determine matrix dimension
D <- max(c(i_vals, j_vals))  # number of antibiotics

# Initialize empty matrix
omega_matrix <- matrix(NA, nrow = D, ncol = D)

# Fill the matrix
for (k in seq_along(omega_hier_means)) {
  i <- i_vals[k]
  j <- j_vals[k]
  omega_matrix[i, j] <- omega_hier_means[k]
  omega_matrix[j, i] <- omega_hier_means[k]  # enforce symmetry
}

antibiotics <- c("Meropenem_bin", 
                 "Piperacillin.tazobactam_bin",
                 "Ampicillin_bin", 
                 "Tigecycline_bin",
                 "Gentamicin_bin", 
                 "Amikacin_bin", 
                 "Ceftazidime_bin",
                 "Levofloxacin_bin", 
                 "Colistin_bin", 
                 "Amoxycillin.clavulanate_bin")  # update as needed
colnames(omega_matrix) <- antibiotics
rownames(omega_matrix) <- antibiotics

ggcorrplot(omega_matrix, 
           type = "lower", 
           lab = TRUE,
           title = "Posterior Mean Correlation Matrix (Omega)")








params_fixed_means <- colMeans(params_fixed)
params_hier_means <- colMeans(params_hier)


X_Y_fixed <- as.data.frame(X_Y_fixed)
X_Y_hier <- as.data.frame(X_Y_hier)
# this results in one point estimate for each parameter. Now we can construct
# multivariate normal distributions for each patient profile

for (param_name in names(params_fixed_means)) {
  X_Y_fixed[[param_name]] <- rep(params_fixed_means[[param_name]], nrow(X_Y_fixed))
}

X_Y_fixed[] <- lapply(X_Y_fixed, function(x) as.numeric(as.character(x)))

covariate_names_fixed <- c("age65 and Over", "genderMale", "icu1",
                     "countryGreece", "countryItaly", "countryPortugal",
                     "countrySlovenia", "countrySpain")


# calculating mu per antibiotic per isolate
K_fixed <- length(covariate_names_fixed)
D <- 10  # Or whatever your number of outcomes is

for (j in 1:D) {
  X_Y_fixed[[paste0("mu_", j)]] <- 
    X_Y_fixed[[paste0("alpha[", j, "]")]] + 
    rowSums(sapply(1:K_fixed, function(k) {
      X_Y_fixed[[covariate_names_fixed[k]]] * X_Y_fixed[[paste0("beta[", k, ",", j, "]")]]
    }))
}







# calculating probability of meropenem resistance conditional on pip/tazo resistance
# for first isolate

antibiotics <- c(
  "Meropenem_bin", 
  "Piperacillin.tazobactam_bin",
  "Ampicillin_bin", 
  "Tigecycline_bin",
  "Gentamicin_bin", 
  "Amikacin_bin", 
  "Ceftazidime_bin",
  "Levofloxacin_bin", 
  "Colistin_bin", 
  "Amoxycillin.clavulanate_bin"
)

n_antibiotics <- length(antibiotics)

# Create a named list to hold a dataframe per conditioning antibiotic
conditional_dfs <- list()

for(cond_i in 1:n_antibiotics){
  
  cond_ab <- antibiotics[cond_i]
  
  # Start with test data (including isolate id)
  df_cond <- X_Y_fixed %>% select(id, all_of(names(X_Y_fixed)[!grepl("^mu_|^Omega", names(X_Y_fixed))]))
  
  for(other_i in 1:n_antibiotics){
    if(other_i != cond_i){
      
      mu_i_col <- paste0("mu_", other_i)
      mu_cond_col <- paste0("mu_", cond_i)
      
      omega_col_name <- if(other_i < cond_i) {
        paste0("Omega[", other_i, ",", cond_i, "]")
      } else {
        paste0("Omega[", cond_i, ",", other_i, "]")
      }
      
      # Check all columns exist
      required_cols <- c(mu_i_col, mu_cond_col, omega_col_name, cond_ab)
      if(!all(required_cols %in% colnames(X_Y_fixed))){
        stop(paste("Missing columns:", paste(setdiff(required_cols, colnames(X_Y_fixed)), collapse=", ")))
      }
      
      # Calculate conditional probs for all rows
      cond_probs <- mapply(function(mu1, mu2, rho, y2_obs){
        mu_vec <- c(mu1, mu2)
        Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
        
        if(y2_obs == 1){
          joint_prob <- pmvnorm(lower = c(0, 0), upper = c(Inf, Inf), mean = mu_vec, sigma = Sigma)[1]
          marginal_prob <- pnorm(mu2)
        } else {
          joint_prob <- pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = mu_vec, sigma = Sigma)[1]
          marginal_prob <- pnorm(-mu2)
        }
        joint_prob / marginal_prob
      },
      mu1 = X_Y_fixed[[mu_i_col]],
      mu2 = X_Y_fixed[[mu_cond_col]],
      rho = X_Y_fixed[[omega_col_name]],
      y2_obs = X_Y_fixed[[cond_ab]])
      
      # Add this conditional probability as a new column
      # Clean antibiotic name for column naming
      other_ab_clean <- gsub("[.]", "_", antibiotics[other_i])
      colname <- paste0("condProb_given_", gsub("[.]", "_", cond_ab), "_for_", other_ab_clean)
      df_cond[[colname]] <- cond_probs
    }
  }
  
  # Store the dataframe in the list, named by conditioning antibiotic
  conditional_dfs[[cond_ab]] <- df_cond
}


# Set output folder
output_folder <- "cond_pred_v2"
dir.create(output_folder, showWarnings = FALSE)

# Loop through each element in the list
for (name in names(conditional_dfs)) {
  # Clean the name and add "_v2" before ".csv"
  clean_name <- gsub("[^A-Za-z0-9_]", "_", name)
  file_path <- file.path(output_folder, paste0("conditional_probs_", clean_name, "_fixed.csv"))
  
  # Save the dataframe
  write.csv(conditional_dfs[[name]], file = file_path, row.names = FALSE)
}









# generating list of all covariate combinations
covariates <- X_south_fixed_test[, -1]  # dropping ID
list_of_combinations <- apply(covariates, 1, function(row) c(1, as.numeric(row)))
list_of_combinations <- split(list_of_combinations, rep(1:ncol(list_of_combinations), each = nrow(list_of_combinations)))
str(list_of_combinations[[1]])  # first isolate's covariates + intercept
head(X_south_fixed_test)

# renaming so that it fits with Nicolas' code
patient_profiles <- list_of_combinations

# vector of antibiotics (preserving order)
antibiotics <- c(
  "Meropenem_bin", 
  "Piperacillin.tazobactam_bin",
  "Ampicillin_bin", 
  "Tigecycline_bin",
  "Gentamicin_bin", 
  "Amikacin_bin", 
  "Ceftazidime_bin",
  "Levofloxacin_bin", 
  "Colistin_bin", 
  "Amoxycillin.clavulanate_bin"
)

get_ab_index <- function(name) {
  match(name, antibiotics)
}

get_conditional_prob_point <- function(i1, i2, beta_mat, omega_mat, x_vec) {
  # i1, i2: outcomes index (1-based)
  # beta_mat: matrix [D, K] (point estimate)
  # omega_mat: matrix [D, D] (point estimate)
  # x_vec: covariate vector (length K)
  
  library(mvtnorm)
  
  mu1 <- sum(beta_mat[i1, ] * x_vec)
  mu2 <- sum(beta_mat[i2, ] * x_vec)
  rho <- omega_mat[i1, i2]
  
  # Joint probability P(Y1=1, Y2=1)
  joint_prob <- pmvnorm(
    lower = c(0, 0),
    upper = c(Inf, Inf),
    mean = c(mu1, mu2),
    sigma = matrix(c(1, rho, rho, 1), nrow = 2)
  )[1]
  
  # Marginal probability P(Y2=1)
  marginal_y2 <- pnorm(mu2)
  
  # Conditional probability P(Y1=1 | Y2=1)
  cond_prob <- joint_prob / marginal_y2
  
  return(cond_prob)
}



results <- list()

for (profile_name in names(patient_profiles)) {
  x_vec <- patient_profiles[[profile_name]]
  
  for (i1 in 1:length(antibiotics)) {
    for (i2 in 1:length(antibiotics)) {
      if (i1 != i2) {
        prob <- get_conditional_prob_point(
          i1 = i1,
          i2 = i2,
          beta_mat = beta_fixed_means,
          omega_mat = omega_fixed_means,
          x_vec = x_vec
        )
        
        results[[length(results) + 1]] <- tibble(
          profile = profile_name,
          ab_i1 = antibiotics[i1],
          ab_i2 = antibiotics[i2],
          cond_prob = prob
        )
      }
    }
  }
}

final_results <- bind_rows(results)

