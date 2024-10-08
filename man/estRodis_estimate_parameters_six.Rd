% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estRodis_estimate_parameters_six.R
\name{estRodis_estimate_parameters_six}
\alias{estRodis_estimate_parameters_six}
\title{Estimate parameters using Model 6}
\usage{
estRodis_estimate_parameters_six(
  clusters_size,
  clusters_freq,
  prior_r = c(10, 10),
  prior_k = c(5, 10),
  prior_mutation = c(27, 68),
  testing_proba = 1,
  sequencing_proba = 1,
  pars = NA,
  chains = 4,
  iter = 2000,
  warmup = floor(iter/2),
  thin = 1,
  seed = sample.int(.Machine$integer.max, 1),
  init = "random",
  check_data = TRUE,
  sample_file = NULL,
  diagnostic_file = NULL,
  verbose = FALSE,
  algorithm = c("NUTS", "HMC", "Fixed_param"),
  control = NULL,
  include = TRUE,
  cores = getOption("mc.cores", 1L),
  open_progress = interactive() && !isatty(stdout()) && !identical(Sys.getenv("RSTUDIO"),
    "1"),
  show_messages = TRUE
)
}
\arguments{
\item{clusters_size}{A list of sizes of identical sequence clusters.}

\item{clusters_freq}{A list of frequencies of identical sequence clusters.}

\item{prior_r}{Shape and rate parameter for the prior distribution (gamma) of
the effective reproduction number.}

\item{prior_k}{Shape and rate parameter for the prior distribution (gamma) of
the dispersion parameter.}

\item{prior_mutation}{Shape parameters for the prior distribution (beta) of the mutation probability.}

\item{testing_proba}{Probability that a case gets confirmed by a test.}

\item{sequencing_proba}{Probability that a case confirmed by a test gets sequenced.}

\item{pars}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{chains}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{iter}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{warmup}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{thin}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{seed}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{init}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{check_data}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{sample_file}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{diagnostic_file}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{verbose}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{algorithm}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{control}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{include}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{cores}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{open_progress}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}

\item{show_messages}{Parameter that is passed on to rstan::sampling, see manual of rstan::sampling for more details.}
}
\value{
An object of S4 class \code{stanfit} containing the fitted results.
}
\description{
\code{estRodis_estimate_parameters_six()} estimates the effective reproduction number and
the dispersion parameter based on the size distribution of identical sequence clusters.
}
\details{
The core of the function \code{estRodis_estimate_parameters_six()} is a mathematical model
of the size distribution of sequence clusters, in which viral transmission,
the mutation of the virus, and incomplete case-detection are integrated.
Parameters are estimated by a Bayesian inference model implemented in Stan.
The effective reproduction number, the dispersion parameter and the mutation probability
are included into the model with prior distributions.
The testing probability and the sequencing probability are included into
the model as fixed parameters.
More details can be found in the following paper: ADD REFERENCE
}
