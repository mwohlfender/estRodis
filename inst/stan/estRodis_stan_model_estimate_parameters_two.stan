
// stan model for estimation of the effective reproduction number, the overdispersion parameter and the number of yearly mutations

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
  // mean time between receving the virus and transmitting it to another person
  real mean_generation_interval;
  // mean and standard deviation for the prior distribution (normal) of the number of yearly mutations
  real prior_number_yearly_mutations[2];
  // probability that a case gets confirmed by a test and gets sequenced
  real detection_proba;
}

parameters {
  real < lower = 0 > R;
  real < lower = 0 > k;
  real < lower = 0 > number_yearly_mutations;
}

transformed parameters {
  real mutation_proba = 1 - exp(- number_yearly_mutations / 365.25 * mean_generation_interval);
}

model {
  // prior distributions
  // alpha (shape) = prior_r[1], beta (rate) = prior_r[2]
  R ~ gamma(prior_r[1], prior_r[2]);
  // alpha (shape) = prior_k[1], beta (rate) = prior_k[2]
  k ~ gamma(prior_k[1], prior_k[2]);
  // mu (mean) = prior_number_yearly_mutations[1], sigma (sd) = prior_number_yearly_mutations[2]
  number_yearly_mutations ~ normal(prior_number_yearly_mutations[1], prior_number_yearly_mutations[2]);

  // likelihood
  target += estRodis_stan_likelihood_log(clusters_size, clusters_freq, R, k, mutation_proba, detection_proba);
}

generated quantities {
  // The posterior predictive distribution
}
