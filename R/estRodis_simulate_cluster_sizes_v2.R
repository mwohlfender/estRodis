
#' @title Simulation of identical sequence clusters.
#'
#' @description `n_clusters` identical sequence clusters of sizes `1` to `max_cluster_size`
#' are simulated. Viral transmission (`R` and `k`), viral mutation
#' (`mutation_proba`) and case detection
#' (`testing_proba` and `sequencing_proba`) are included in the simulation.
#'
#' @param n_clusters Number of clusters to be simulated.
#' @param max_cluster_size Maximal size of a simulated identical sequence cluster.
#' @param R Effective reproduction number.
#' @param k Dispersion parameter.
#' @param mutation_proba Probability that a mutation occurs at a case.
#' @param testing_proba Probability that a case is confirmed by a test.
#' @param sequencing_proba Probability that the genome of a confirmed case is sequenced.
#'
#' @details More details can be found in the following paper: ADD REFERENCE
#'
#' @return A `data.frame` with three columns, `size`, `frequency` and `percentage`,
#' containing size and frequency of the simulated identical sequence clusters as well as the
#' proportion of the identical sequence clusters of each size.
#'
#' @export
#'
#' @examples estRodis_simulate_cluster_sizes()

estRodis_simulate_cluster_sizes_v2 <- function(n_clusters = 1000,
                                               max_cluster_size = 2500,
                                               R = 1.0,
                                               k = 0.3,
                                               mutation_proba = 0.2,
                                               testing_proba = 0.6,
                                               sequencing_proba = 0.4) {

  # calculate detection probability
  detection_proba <- testing_proba * sequencing_proba

  # create data frame to store the size of the simulated identical sequence clusters
  n_nodes_identical_sequence_clusters_sim <- data.frame(size = c(0:max_cluster_size),
                                                        frequency = rep(x = 0, times = max_cluster_size + 1))

  # simulate identical sequence clusters until at least n_clusters identical sequence clusters have been created
  # convention: the size of an identical sequence cluster is the number of detected nodes that have transmitted the respective variant to their (direct) offspring (might have zero offspring)
  while (sum(n_nodes_identical_sequence_clusters_sim |> dplyr::filter(.data$size >= 1) |> dplyr::pull(.data$frequency)) < n_clusters) {

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

  return(n_nodes_identical_sequence_clusters_sim |> dplyr::filter(.data$frequency != 0, .data$size >= 1) |> dplyr::mutate(percentage = .data$frequency / sum(.data$frequency)))

}
