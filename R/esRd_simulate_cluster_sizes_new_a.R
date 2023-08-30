
#' @title Simulation of identical sequence clusters.
#'
#' @description `n_clusters` identical sequence clusters of sizes `1` to `max_cluster_size`
#' are simulated. Viral transmission (`R` and `k`), viral mutation
#' (`yearly_mutation_rate` and `mean_generation_interval`) and case detection
#' (`testing_proba` and `sequencing_proba`) are included in the simulation.
#'
#' @param n_clusters Number of clusters to be simulated.
#' @param max_cluster_size Maximal size of a simulated identical sequence cluster.
#' @param R Effective reproduction number.
#' @param k Dispersion parameter.
#' @param yearly_mutation_rate Number of mutations of the viral genome per year.
#' @param mean_generation_interval Mean generation interval.
#' @param testing_proba Probability that a case is confirmed by a test.
#' @param sequencing_proba Probability that the genome of a confirmed case is sequenced.
#'
#' @details
#'
#' @return
#'
#' @export
#'
#' @examples esRd_simulate_cluster_sizes()
#'
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom stats rbinom
#' @importFrom stats rnbinom
esRd_simulate_cluster_sizes_new_a <- function(n_clusters = 1000,
                                              max_cluster_size = 2500,
                                              R = 1.0,
                                              k = 0.3,
                                              yearly_mutation_rate = 14,
                                              mean_generation_interval = 5.2,
                                              testing_proba = 0.6,
                                              sequencing_proba = 0.4) {

  # calculate mutation probability
  mutation_proba <- 1 - exp(- yearly_mutation_rate / 365.25 * mean_generation_interval)

  # calculate detection probability
  detection_proba <- testing_proba * sequencing_proba

  # create data frame to store the size of the simulated identical sequence clusters
  n_nodes_identical_sequence_clusters_sim <- data.frame(size = c(0:max_cluster_size),
                                                        frequency = rep(x = 0, times = max_cluster_size + 1))

  # simulate identical sequence clusters until at least n_clusters identical sequence clusters have been created
  # convention: the size of an identical sequence cluster is the number of detected nodes that have transmitted the respective variant to their (direct) offspring (might have zero offspring)
  while (sum(n_nodes_identical_sequence_clusters_sim |> dplyr::filter(size >= 1) |> dplyr::pull(frequency)) < n_clusters) {

    # create initial case and determine whether it is detected or not
    n_free_leaves <- 1
    n_nodes_detected <- stats::rbinom(n = 1, size = 1, prob = detection_proba)

    # add new cases to the identical sequence cluster as long as there are free leaves,
    # i.e. cases for which the number of offspring has not been determined yet
    while ((n_free_leaves > 0) && (n_nodes_detected <= max_cluster_size)) {

      # number of offspring of the free leaf
      n_new_leaves_tree <- stats::rnbinom(n = 1, size = k, mu = R)

      # number of offspring in the same identical sequence cluster of the free leaf
      n_new_leaves_cluster <- sum(stats::rbinom(n = n_new_leaves_tree, size = 1, prob = 1 - mutation_proba))

      # determine which of the new cases are detected and update the number of detected cases of the identical sequence cluster
      n_nodes_detected <- n_nodes_detected + sum(stats::rbinom(n = n_new_leaves_cluster, size = 1, prob = detection_proba))

      # update the number of free leaves of the identical sequence cluster
      # add the new cases of the identical sequence cluster and remove the free leaf whose offspring have been determined
      n_free_leaves <- n_free_leaves + n_new_leaves_cluster - 1

    }

    # store the size of the identical sequence cluster
    n_nodes_identical_sequence_clusters_sim$frequency[min(n_nodes_detected, max_cluster_size) + 1] <- n_nodes_identical_sequence_clusters_sim$frequency[min(n_nodes_detected, max_cluster_size) + 1] + 1

  }

  return(n_nodes_identical_sequence_clusters_sim |> dplyr::filter(frequency != 0, size >= 1) |> dplyr::mutate(percentage = frequency / sum(frequency)))

}


