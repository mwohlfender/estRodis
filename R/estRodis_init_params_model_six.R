
#' @title Initial values for Model 6
#'
#' @description Sample from prior distributions of R, k and the mutation probability.
#'
#' @param prior_r Parameters for prior distribution of R (shape and rate parameter for gamma distribution)
#' @param prior_k Parameters for prior distribution of k (shape and rate parameter for gamma distribution)
#' @param prior_mutation Parameters for prior distributions of mutation probability (shape parameters for beta distribution)
#'
#' @return A named list containing samples of the prior distributions of R, k and the mutation probability.
#'
#' @export
#'
#' @examples
#' estRodis_init_params_model_six()

estRodis_init_params_model_six <- function(prior_r = c(10, 10), prior_k = c(5, 10), prior_mutation = c(27, 68)) {

  # initial value for R sampled from gamma distribution(prior_r[1], prior_r[2])
  # initial value for k sampled from gamma distribution(prior_k[1], prior_k[2])
  # initial value for mutation_proba sampled from beta distribution(prior_mutation[1], prior_mutation[2])
  result <- list(R = stats::rgamma(1, prior_r[1], prior_r[2]),
                 k = stats::rgamma(1, prior_k[1], prior_k[2]),
                 mutation_proba = stats::rbeta(n = 1, shape1 = prior_mutation[1], shape2 = prior_mutation[2]))

  return(result)

}
