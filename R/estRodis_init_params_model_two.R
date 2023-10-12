
#' @title Initial values for Model 2
#'
#' @description Sample from prior distributions of R, k, the number of yearly mutations and the testing probability.
#'
#' @param prior_r Parameters for prior distribution of R (gamma)
#' @param prior_k Parameters for prior distribution of k (gamma)
#' @param prior_number_yearly_mutations Parameters for prior distribution of the number of yearly mutations (normal)
#'
#' @return A named list containing samples of the prior distributions of R, k and the number of yearly mutations.
#'
#' @export
#'
#' @examples
#' estRodis_init_params_model_two()

estRodis_init_params_model_two <- function(prior_r = c(10, 10), prior_k = c(5, 10), prior_number_yearly_mutations = c(14, 0.5)) {

  # initial value for R sampled from gamma distribution(prior_r[1], prior_r[2])
  # initial value for k sampled from gamma distribution(prior_k[1], prior_k[2])
  # initial value for number_yearly_mutations sampled from normal distribution(prior_number_yearly_mutations[1], prior_number_yearly_mutations[2])
  result <- list(R = stats::rgamma(1, prior_r[1], prior_r[2]),
                 k = stats::rgamma(1, prior_k[1], prior_k[2]),
                 number_yearly_mutations = stats::rnorm(1, prior_number_yearly_mutations[1], prior_number_yearly_mutations[2]))

  return(result)

}
