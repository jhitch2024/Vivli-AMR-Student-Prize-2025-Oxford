data {
  int<lower=1> N;                            // Number of isolates
  int<lower=1> P;                            // Number of covariates
  int<lower=1> D;                            // Number of outcomes (e.g., antibiotics)
  matrix[N, P] X;                            // Covariate matrix
  array[N, D] int<lower=0, upper=1> Y;       // Binary resistance outcomes
  int<lower=1> C;                            // Number of countries
  array[N] int<lower=1, upper=C> country_id; // Country ID for each isolate

  // Test data
  int<lower=1> N_test;                       // Number of test isolates
  matrix[N_test, P] X_test;                  // Test covariate matrix
  array[N_test] int<lower=1, upper=C> country_id_test; // Country ID for test isolates
}

parameters {
  matrix[C, D] alpha_raw;                    // Raw country intercepts
  vector[D] mu_alpha;                        // Mean for country intercepts
  vector<lower=0>[D] sigma_alpha;            // Std dev for country intercepts
  matrix[P, D] beta;                         // Covariate effects
  cholesky_factor_corr[D] L;                 // Cholesky factor of correlation matrix
  matrix[N, D] Z_raw;                        // Raw latent variables for training
}

transformed parameters {
  matrix[C, D] alpha;                        // Country-level intercepts
  matrix[N, D] Z;                            // Latent variables for training

  for (c in 1:C)
    for (d in 1:D)
      alpha[c, d] = mu_alpha[d] + sigma_alpha[d] * alpha_raw[c, d];

  for (n in 1:N)
    Z[n] = alpha[country_id[n]] + X[n] * beta + Z_raw[n] * L';
}

model {
  // Priors
  to_vector(alpha_raw) ~ normal(0, 1);
  mu_alpha ~ normal(0, 1);
  sigma_alpha ~ exponential(1);
  to_vector(beta) ~ normal(0, 1);
  L ~ lkj_corr_cholesky(4);
  to_vector(Z_raw) ~ normal(0, 1);

  // Likelihood
  for (n in 1:N)
    for (d in 1:D)
      Y[n, d] ~ bernoulli(Phi(Z[n, d]));
}

generated quantities {
  array[N, D] int Y_rep;                    // Posterior predictive replicates for training
  matrix[N, D] log_lik;                     // Log-likelihood for training
  matrix[N, D] pred_prob;                   // Predicted probabilities for training
  corr_matrix[D] Omega;                     // Correlation matrix

  matrix[N_test, D] Z_test_raw;             // Raw latent vars for test data
  matrix[N_test, D] Z_test;                 // Latent vars for test data
  matrix[N_test, D] pred_prob_test;         // Predicted probabilities test data
  array[N_test, D] int Y_test_rep;          // Predictive replicates test data

  Omega = multiply_lower_tri_self_transpose(L);

  // Training predictions and log likelihood
  for (n in 1:N) {
    for (d in 1:D) {
      real p = Phi(Z[n, d]);
      pred_prob[n, d] = p;
      Y_rep[n, d] = bernoulli_rng(p);
      log_lik[n, d] = bernoulli_lpmf(Y[n, d] | p);
    }
  }

  // Test predictions incorporating correlations and country intercepts
  for (n in 1:N_test) {
    for (d in 1:D) {
      Z_test_raw[n, d] = normal_rng(0, 1);
    }
    Z_test[n] = alpha[country_id_test[n]] + X_test[n] * beta + Z_test_raw[n] * L';

    for (d in 1:D) {
      real p = Phi(Z_test[n, d]);
      pred_prob_test[n, d] = p;
      Y_test_rep[n, d] = bernoulli_rng(p);
    }
  }
}
