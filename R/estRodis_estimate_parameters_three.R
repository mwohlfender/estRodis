
#' @title Estimate parameters using Model 3
#'
#' @description `estRodis_estimate_parameters_three()` estimates the effective reproduction number,
#' the dispersion parameter and the probability that a case is confirmed by a test
#' based on the size distribution of identical sequence clusters.
#'
#' @param clusters_size A list of sizes of identical sequence clusters.
#' @param clusters_freq A list of frequencies of identical sequence clusters.
#' @param prior_r Shape and rate parameter for the prior distribution (gamma) of
#' the effective reproduction number.
#' @param prior_k Shape and rate parameter for the prior distribution (gamma) of
#' the dispersion parameter.
#' @param mutation_proba Probability that a case undergoes a mutation.
#' @param prior_testing Parameters (shape and interval boundaries) for the prior
#' distribution (scaled beta) of the fraction of cases that are confirmed by a test.
#' @param sequencing_proba Probability that a case confirmed by a test gets sequenced.
#' @param pars Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param chains Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param iter Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param warmup Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param thin Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param seed Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param init Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param check_data Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param sample_file Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param diagnostic_file Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param verbose Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param algorithm Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param control Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param include Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param cores Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param open_progress Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#' @param show_messages Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.
#'
#' @details The core of the function `estRodis_estimate_parameters_three()` is a mathematical model
#' of the size distribution of sequence clusters, in which viral transmission,
#' the mutation of the virus, and incomplete case-detection are integrated.
#' Parameters are estimated by a Bayesian inference model implemented in Stan.
#' The effective reproduction number, the dispersion parameter, the number of yearly mutations
#' and the testing probability are included into the model with prior distributions.
#' The mean generation interval and the sequencing probability are included into
#' the model as fixed parameters.
#' More details can be found in the following paper: ADD REFERENCE
#'
#' @return An object of S4 class `stanfit` containing the fitted results.
#'
#' @export


estRodis_estimate_parameters_three <- function(clusters_size,
                                               clusters_freq,
                                               prior_r = c(10, 10),
                                               prior_k = c(5, 10),
                                               mutation_proba = 0.2811,
                                               prior_testing = c(1, 3, 0.05, 1),
                                               sequencing_proba = 1,
                                               pars = NA,
                                               chains = 4,
                                               iter = 2000,
                                               warmup = floor(iter/2),
                                               thin = 1,
                                               seed = sample.int(.Machine$integer.max, 1),
                                               init = 'random',
                                               check_data = TRUE,
                                               sample_file = NULL,
                                               diagnostic_file = NULL,
                                               verbose = FALSE,
                                               algorithm = c("NUTS", "HMC", "Fixed_param"),
                                               control = NULL,
                                               include = TRUE,
                                               cores = getOption("mc.cores", 1L),
                                               open_progress = interactive() && !isatty(stdout()) && !identical(Sys.getenv("RSTUDIO"), "1"),
                                               show_messages = TRUE) {

  standata <- list(M = length(clusters_size),
                   clusters_size = clusters_size,
                   clusters_freq = clusters_freq,
                   prior_r = prior_r,
                   prior_k = prior_k,
                   mutation_proba = mutation_proba,
                   prior_testing = prior_testing,
                   sequencing_proba = sequencing_proba)

  out <- rstan::sampling(object = stanmodels$estRodis_stan_model_estimate_parameters_three,
                         data = standata,
                         pars = pars,
                         chains = chains,
                         iter = iter,
                         warmup = warmup,
                         thin = thin,
                         seed = seed,
                         init = init,
                         check_data = check_data,
                         sample_file = sample_file,
                         diagnostic_file = diagnostic_file,
                         verbose = verbose,
                         algorithm = algorithm,
                         control = control,
                         include = include,
                         cores = cores,
                         open_progress = open_progress,
                         show_messages = show_messages)

  return(out)

}


