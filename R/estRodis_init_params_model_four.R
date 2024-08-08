
#' @title Initial values for Model 4
#'
#' @description Sample from prior distributions of R and k.
#'
#' @param prior_r Parameters for prior distribution of R (shape and rate parameter for gamma distribution)
#' @param prior_k Parameters for prior distribution of k (shape and rate parameter for gamma distribution)
#'
#' @return A named list containing samples of the prior distributions of R and k.
#'
#' @export
#'
#' @examples
#' estRodis_init_params_model_four()

estRodis_init_params_model_four <- function(prior_r = c(10, 10), prior_k = c(5, 10)) {

  # initial value for R sampled from gamma distribution(shape = prior_r[1], rate = prior_r[2])
  # initial value for k sampled from gamma distribution(shape = prior_k[1], rate = prior_k[2])
  result <- list(R = stats::rgamma(n = 1, shape = prior_r[1], rate = prior_r[2]),
                 k = stats::rgamma(n = 1, shape = prior_k[1], rate = prior_k[2]))

  return(result)

}
