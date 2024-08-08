
#' @title Initial values for Model 1
#'
#' @description Sample from prior distributions of R, k, the number of yearly mutations and the testing probability.
#'
#' @param prior_r Parameters for prior distribution of R (shape and rate parameter for gamma distribution)
#' @param prior_k Parameters for prior distribution of k (shape and rate parameter for gamma distribution)
#' @param prior_number_yearly_mutations Parameters for prior distribution of the number of yearly mutations (mean and standard deviation for normal distribution)
#' @param prior_testing Parameters for prior distribution of testing probability (scaled beta distribution)
#'
#' @return A named list containing samples of the prior distributions of R, k, the number of yearly mutations and the testing probability.
#'
#' @export
#'
#' @examples
#' estRodis_init_params_model_one()

estRodis_init_params_model_one <- function(prior_r = c(10, 10), prior_k = c(5, 10), prior_number_yearly_mutations = c(14, 0.5), prior_testing = c(1, 3, 0.05, 1)) {

  distribution <- unlist(lapply(X = seq(prior_testing[3], prior_testing[4], 0.00001),
                                FUN = function(x) estRodis_scaled_beta_distribution_pdf(prior_testing[1], prior_testing[2], prior_testing[3], prior_testing[4], x)))

  # initial value for R sampled from gamma distribution(shape = prior_r[1], rate = prior_r[2])
  # initial value for k sampled from gamma distribution(shape = prior_k[1], rate = prior_k[2])
  # initial value for number_yearly_mutations sampled from normal distribution(mean = prior_number_yearly_mutations[1], sd = prior_number_yearly_mutations[2])
  # initial value for testing_proba sampled from scaled beta distribution(prior_testing[1], prior_testing[2]) on the interval (prior_testing[3], prior_testing[4])
  result <- list(R = stats::rgamma(n = 1, shape = prior_r[1], rate = prior_r[2]),
                 k = stats::rgamma(n = 1, shape = prior_k[1], rate = prior_k[2]),
                 number_yearly_mutations = stats::rnorm(n = 1, mean = prior_number_yearly_mutations[1], sd = prior_number_yearly_mutations[2]),
                 testing_proba =  sample(seq(prior_testing[3], prior_testing[4], 0.00001), 1, replace = TRUE, prob = distribution))

  return(result)

}
