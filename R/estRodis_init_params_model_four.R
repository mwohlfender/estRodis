
#' @title Initial values for Model 4
#'
#' @description Sample from prior distributions of R and k.
#'
#' @param prior_r Parameters for prior distribution of R (gamma)
#' @param prior_k Parameters for prior distribution of k (gamma)
#'
#' @return A named list containing samples of the prior distributions of R and k.
#'
#' @export
#'
#' @examples
#' estRodis_init_params_model_four()

estRodis_init_params_model_four <- function(prior_r = c(10, 10), prior_k = c(5, 10)) {

  # initial value for R sampled from gamma distribution(prior_r[1], prior_r[2])
  # initial value for k sampled from gamma distribution(prior_k[1], prior_k[2])
  result <- list(R = stats::rgamma(1, prior_r[1], prior_r[2]),
                 k = stats::rgamma(1, prior_k[1], prior_k[2]))

  return(result)

}
