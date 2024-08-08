
#' @title Initial values for Model 5
#'
#' @description Sample from prior distributions of R, k, the mutation probability and the testing probability.
#'
#' @param prior_r Parameters for prior distribution of R (shape and rate parameter for gamma distribution)
#' @param prior_k Parameters for prior distribution of k (shape and rate parameter for gamma distribution)
#' @param prior_mutation Parameters for prior distributions of mutation probability (shape parameters for beta distribution)
#' @param prior_testing Parameters for prior distribution of testing probability (scaled beta distribution)
#'
#' @return A named list containing samples of the prior distributions of R, k, the mutation probability and the testing probability.
#'
#' @export
#'
#' @examples
#' estRodis_init_params_model_five()

estRodis_init_params_model_five <- function(prior_r = c(10, 10), prior_k = c(5, 10), prior_mutation = c(27, 68), prior_testing = c(1, 3, 0.05, 1)) {

  distribution <- unlist(lapply(X = seq(prior_testing[3], prior_testing[4], 0.00001),
                                FUN = function(x) estRodis_scaled_beta_distribution_pdf(prior_testing[1], prior_testing[2], prior_testing[3], prior_testing[4], x)))

  # initial value for R sampled from gamma distribution(prior_r[1], prior_r[2])
  # initial value for k sampled from gamma distribution(prior_k[1], prior_k[2])
  # initial value for mutation_proba sampled from beta distribution(prior_mutation[1], prior_mutation[2])
  # initial value for testing_proba sampled from scaled beta distribution(prior_testing[1], prior_testing[2]) on the interval (prior_testing[3], prior_testing[4])
  result <- list(R = stats::rgamma(n = 1, prior_r[1], prior_r[2]),
                 k = stats::rgamma(n = 1, prior_k[1], prior_k[2]),
                 mutation_proba = stats::rbeta(n = 1, shape1 = prior_mutation[1], shape2 = prior_mutation[2]),
                 testing_proba =  sample(seq(prior_testing[3], prior_testing[4], 0.00001), 1, replace = TRUE, prob = distribution))

  return(result)

}
