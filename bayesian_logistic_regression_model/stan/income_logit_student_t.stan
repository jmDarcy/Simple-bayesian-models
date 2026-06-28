data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  int<lower=0, upper=1> y[N];
  vector[K] beta_prior_location;
  vector<lower=0>[K] beta_prior_scale;
  real<lower=1> beta_prior_df;
  int<lower=1> N_pred;
  matrix[N_pred, K] X_pred;
}

parameters {
  vector[K] beta;
}

model {
  beta ~ student_t(beta_prior_df, beta_prior_location, beta_prior_scale);
  y ~ bernoulli_logit(X * beta);
}

generated quantities {
  vector[N] eta;
  vector[N] p_hat;
  vector[N] log_lik;
  int<lower=0, upper=1> y_rep[N];
  vector[N_pred] p_pred;
  int<lower=0, upper=1> y_pred_rep[N_pred];

  eta = X * beta;

  for (n in 1:N) {
    p_hat[n] = inv_logit(eta[n]);
    log_lik[n] = bernoulli_logit_lpmf(y[n] | eta[n]);
    y_rep[n] = bernoulli_rng(p_hat[n]);
  }

  for (m in 1:N_pred) {
    p_pred[m] = inv_logit(dot_product(row(X_pred, m), beta));
    y_pred_rep[m] = bernoulli_rng(p_pred[m]);
  }
}
