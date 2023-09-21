

#' @title Probability density function of scaled Beta distribution
#'
#' @description `estRodis_scaled_beta_distribution_pdf` evaluates the probability density function of the scaled Beta distribution that is
#' determined by the parameters `alpha` and `beta` and defined on the interval  `p`, `q` at `x`.
#'
#' @param alpha shape parameter (real number > 0)
#' @param beta shape parameter (real number > 0)
#' @param p lower boundary of interval on which the distribution is defined
#' @param q upper boundary of interval on which the distribution is defined
#' @param x list of values at which the probability density function of the scaled beta distribution shall be evaluated
#'
#' @details Values of `x` outside the interval `p`, `q` are evaluated to `0`.
#'
#' @return A list of real numbers.
#' @export
#'
estRodis_scaled_beta_distribution_pdf <- function(alpha, beta, p, q, x) {

  result <- unlist(lapply(X = x,
                          FUN = function(x)
                            if (x >= p & x <= q)
                              {((x-p)^(alpha-1) * (q-x)^(beta-1)) / ((q-p)^(alpha+beta-1) * gamma(alpha) * gamma(beta) / gamma(alpha+beta))}
                            else
                              { 0 }))

  return(result)

}
