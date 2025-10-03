data {
  int<lower=1> N;                            // Number of training samples
  int<lower=1> P;                            // Number of covariates
  int<lower=1> D;                            // Number of outcomes (e.g., antibiotics)
  matrix[N, P] X;                            // Covariate matrix
  array[N, D] int<lower=0, upper=1> Y;       // Binary resistance outcomes

  // Test data
  int<lower=1> N_test;
  matrix[N_test, P] X_test;
}

parameters {
  vector[D] alpha;                           // Intercepts
  matrix[P, D] beta;                         // Covariate effects
  cholesky_factor_corr[D] L;                // Cholesky factor of correlation matrix
  matrix[N, D] Z_raw;                        // Raw latent variables for training
}

transformed parameters {
  matrix[N, D] Z;

  for (n in 1:N)
    Z[n] = alpha' + X[n] * beta + Z_raw[n] * L';
}

model {
  // Priors
  alpha ~ normal(0, 1);
  to_vector(beta) ~ normal(0, 1);
  L ~ lkj_corr_cholesky(4);
  to_vector(Z_raw) ~ normal(0, 1);

  // Likelihood
  for (n in 1:N)
    for (d in 1:D)
      Y[n, d] ~ bernoulli(Phi(Z[n, d]));
}

generated quantities {
  array[N, D] int Y_rep;
  matrix[N, D] log_lik;
  matrix[N, D] pred_prob;
  corr_matrix[D] Omega;

  // Test predictions
  array[N_test, D] int Y_test_rep;
  matrix[N_test, D] pred_prob_test;
  matrix[N_test, D] Z_test_raw;
  matrix[N_test, D] Z_test;

  Omega = multiply_lower_tri_self_transpose(L);

  // Training predictions
  for (n in 1:N) {
    for (d in 1:D) {
      real p = Phi(Z[n, d]);
      pred_prob[n, d] = p;
      Y_rep[n, d] = bernoulli_rng(p);
      log_lik[n, d] = bernoulli_lpmf(Y[n, d] | p);
    }
  }

  // Test latent variables and predictions (incorporating correlations!)
for (n in 1:N_test) {
  for (d in 1:D) {
    Z_test_raw[n, d] = normal_rng(0, 1);
  }
}

Z_test = rep_matrix(alpha', N_test) + X_test * beta + Z_test_raw * L';

for (n in 1:N_test) {
  for (d in 1:D) {
    real p = Phi(Z_test[n, d]);
    pred_prob_test[n, d] = p;
    Y_test_rep[n, d] = bernoulli_rng(p);
  }
}
}
