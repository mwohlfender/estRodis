
#' @title Estimate parameters.
#'
#' @description `esRd_estimate_parameters()` estimates the effective reproduction number,
#' the dispersion parameter, the probability of a case undergoing a mutation and
#' the probability that a case is confirmed by a test based on the size distribution
#' of identical sequence clusters.
#'
#' @param clusters_size A list of sizes of identical sequence clusters.
#' @param clusters_freq A list of frequencies of identical sequence clusters.
#' @param prior_r Shape and rate parameter for the prior distribution (gamma) of
#' the effective reproduction number.
#' @param prior_k Shape and rate parameter for the prior distribution (gamma) of
#' the dispersion parameter.
#' @param mean_generation_interval Mean generation interval.
#' @param prior_number_yearly_mutations Mean and standard deviation for the prior
#' distribution (normal) of the number of yearly mutations.
#' @param prior_testing Parameters (shape and interval boundaries) for the prior
#' distribution (scaled beta) of the fraction of cases that are confirmed by a test.
#' @param sequencing_proba Sequencing probability.
#'
#' @details The core of the function `esRd_estimate_parameters()` is a mathematical model
#' of the size distribution of sequence clusters, in which viral transmission,
#' the mutation of the virus, and incomplete case-detection are integrated.
#' Parameters are estimated by a Bayesian inference model implemented in Stan.
#'
#' @return An object of S4 class `stanfit` containing the fitted results.
#'
#' @export
#'
#' @examples
#' options(mc.cores = parallel::detectCores())
#' esRd_estimate_parameters(
#'   c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 20, 22, 23, 28, 29, 32, 33, 34, 35, 36, 83, 103),
#'   c(703, 117, 49, 37, 19, 17, 5, 15, 4, 3, 1, 3, 4, 2, 2, 1, 3, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1),
#'   c(10,10),
#'   c(5,10),
#'   5.2,
#'   c(14,0.5),
#'   c(1,3,0.05,1),
#'   1)

esRd_estimate_parameters <- function(clusters_size,
                                     clusters_freq,
                                     prior_r = c(10, 10),
                                     prior_k = c(5, 10),
                                     mean_generation_interval  = 5.2,
                                     prior_number_yearly_mutations = c(14, 0.5),
                                     prior_testing = c(1, 3, 0.05, 1),
                                     sequencing_proba,
                                     chains = 4,
                                     iter = 2000,
                                     warmup = floor(iter/2),
                                     thin = 1,
                                     control = NULL) {

  standata <- list(M = length(clusters_size),
                   clusters_size = clusters_size,
                   clusters_freq = clusters_freq,
                   prior_r = prior_r,
                   prior_k = prior_k,
                   mean_generation_interval = mean_generation_interval,
                   prior_number_yearly_mutations = prior_number_yearly_mutations,
                   prior_testing = prior_testing,
                   sequencing_proba = sequencing_proba)

  out <- rstan::sampling(object = stanmodels$stan_model_R_k_testing_v1_a,
                         data = standata,
                         chains = chains, iter = iter, warmup = warmup, thin = thin, control = control)

  return(out)

}

