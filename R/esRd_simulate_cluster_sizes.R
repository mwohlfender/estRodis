

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
#' @importFrom dplyr pull
esRd_simulate_cluster_sizes <- function(n_clusters = 1000,
                                        max_cluster_size = 2500,
                                        R = 1.0,
                                        k = 0.3,
                                        yearly_mutation_rate = 14,
                                        mean_generation_interval = 5.2,
                                        testing_proba = 0.6,
                                        sequencing_proba = 0.4) {

  # distribution of the size of identical sequence clusters
  n_nodes_identical_sequence_clusters_sim <- data.frame(size = c(0:max_cluster_size),
                                                        frequency = rep(x = 0, times = max_cluster_size + 1))

  # simulate transmission trees and store the sizes of the identical sequence clusters they contain
  # until at least n_clusters identical sequence clusters have been created
  # convention: the size of an identical sequence cluster is the number of nodes that have transmitted the respective variant to their (direct) offspring (might have zero offspring)
  while (sum(n_nodes_identical_sequence_clusters_sim$frequency[-1]) < n_clusters) {

    # 1 - simulate transmission tree ----
    # offspring distribution: negative binomial distribution with parameters size = k (dispersion parameter) and mu = R (mean)
    tree <- esRd_create_tree_structure(R = R, k = k, max_tree_size = max_cluster_size)

    nodes <- tree[[1]]
    leaves <- tree[[2]]
    edges <- tree[[3]]

    # 2 - simulate mutations ----
    # whether a mutation occurs at a node or not is determined by a Bernoulli distribution with parameter mutation_proba
    mutation_proba <- 1 - exp(- yearly_mutation_rate / 365.25 * mean_generation_interval)

    tree_mutation_information <- esRd_apply_mutations(nodes = nodes, edges = edges, mutation_proba = mutation_proba)

    nodes_a <- tree_mutation_information[[1]]
    edges_a <- tree_mutation_information[[2]]

    # 3 - simulate detection ----
    # whether a node is detected or not is determined by a Bernoulli distribution with parameter detection_proba
    detection_proba <- testing_proba * sequencing_proba

    tree_detection <- esRd_apply_case_detection(nodes = nodes_a, edges = edges_a, detection_proba = detection_proba)

    nodes_b <- tree_detection[[1]]
    edges_b <- tree_detection[[2]]

    # determine the variants present in the tree whose corresponding identical sequence cluster have size at least 1:
    # variant 0, the variant that the first node received, is such a variant if no mutation occurred at the first node
    # remark: the other variants (all variant except variant 0) are named as follows: number of node at which they appeared as the result of a mutation
    if (nodes_a$mutation_occured[1] == 1) {

      variants <- nodes_a |> dplyr::filter(mutation_occured == 1) |> dplyr::pull(node_key)

    } else {

      variants <- c(0, nodes_a |> dplyr::filter(mutation_occured == 1) |> dplyr::pull(node_key))

    }

    # go through all variants
    for (jj in variants) {

      # check whether some clusters still need to be simulated
      if (sum(n_nodes_identical_sequence_clusters_sim$frequency[-1]) < n_clusters) {

        # determine nodes belonging to the identical sequence cluster
        variant_subtree_nodes <- nodes_a |> dplyr::filter((mutation_occured == 0 & variant_received == jj) | (mutation_occured == 1 & variant_after_mutation == jj))

        # check whether one (or more) of the nodes of the variant subtree is a free leaf
        if (length(intersect(variant_subtree_nodes$node_key, leaves$node_key)) == 0) {

          # determine the size of the identical sequence cluster (random detection with probability detection_proba is applied) corresponding to the variant we are currently looking at
          n_nodes_detected <- nrow(nodes_b |> dplyr::filter((mutation_occured == 0 & variant_received == jj) | (mutation_occured == 1 & variant_after_mutation == jj)) |> dplyr::filter(detection == 1))

          # store the size of the identical sequence cluster (random detection with probability detection_proba is applied) corresponding to the variant we are currently looking at
          n_nodes_identical_sequence_clusters_sim$frequency[min(n_nodes_detected, max_cluster_size) + 1] <-  n_nodes_identical_sequence_clusters_sim$frequency[min(n_nodes_detected, max_cluster_size) + 1] + 1

        } else {

          # complete the simulation of a variant subtree
          completed_variant_subtree <- esRd_complete_variant_subtree(tree_nodes = nodes_b, tree_leaves = leaves, tree_edges = edges_b, variant = jj, limit_size = max_cluster_size, R = R, k = k, mutation_proba = mutation_proba, detection_proba = detection_proba)

          # determine the size of the identical sequence cluster (random detection with probability detection_proba is applied) corresponding to the variant we are currently looking at
          n_nodes_detected <- nrow(completed_variant_subtree[[1]] |> dplyr::filter(detection == 1))

          # store the size of the identical sequence cluster (random detection with probability detection_proba is applied) corresponding to the variant we are currently looking at
          n_nodes_identical_sequence_clusters_sim$frequency[min(n_nodes_detected, max_cluster_size) + 1] <- n_nodes_identical_sequence_clusters_sim$frequency[min(n_nodes_detected, max_cluster_size) + 1] + 1

        }

      }

    }

  }

  return(n_nodes_identical_sequence_clusters_sim |> dplyr::filter(frequency != 0))

}
