
// stan model for estimation of the effective reproduction number, the overdispersion parameter and the testing probability

#include "/stan_functions/estRodis_stan_functions_estimate_parameters.stan"

data {
  // number of different identical sequence cluster sizes
  int < lower = 1 > M;
  // sizes of identical sequence clusters
  int clusters_size[M];
  // frequencies of identical sequence cluster sizes
  int clusters_freq[M];
  // shape and rate parameter for the prior distribution (gamma) of the effective reproduction number
  real prior_r[2];
  // shape and rate parameter for the prior distribution (gamma) of the dispersion parameter
  real prior_k[2];
  // mutation probability
  real mutation_proba;
  // parameters (shape and interval boundaries) for the prior distribution (scaled beta) of
  // the fraction of cases that are confirmed by a test
  real prior_testing[4];
  // probability that a case confirmed by a test gets sequenced
  real sequencing_proba;
}

parameters {
  real < lower = 0 > R;
  real < lower = 0 > k;
  real < lower = prior_testing[3], upper = prior_testing[4]> testing_proba;
}

transformed parameters {
  real detection_proba = testing_proba * sequencing_proba;
}

model {
  // prior distributions
  // alpha (shape) = prior_r[1], beta (rate) = prior_r[2]
  R ~ gamma(prior_r[1], prior_r[2]);
  // alpha (shape) = prior_k[1], beta (rate) = prior_k[2]
  k ~ gamma(prior_k[1], prior_k[2]);
  // alpha = prior_testing[1], beta = prior_testing[2], p = prior_testing[3], q = prior_testing[4]
  testing_proba ~ estRodis_stan_scaled_beta(prior_testing[1], prior_testing[2], prior_testing[3], prior_testing[4]);

  // likelihood
  target += estRodis_stan_likelihood_log(clusters_size, clusters_freq, R, k, mutation_proba, detection_proba);
}

generated quantities {
  // The posterior predictive distribution
}
