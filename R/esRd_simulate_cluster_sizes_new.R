
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom stats rbinom
#' @importFrom stats rnbinom
esRd_simulate_cluster_sizes_new <- function(n_clusters = 1000,
                                            max_cluster_size = 2500,
                                            R = 1.0,
                                            k = 0.3,
                                            yearly_mutation_rate = 14,
                                            mean_generation_interval = 5.2,
                                            testing_proba = 0.6,
                                            sequencing_proba = 0.4) {

  # create data frame to store information about the nodes of the variant subtree we consider
  # we define it to be most likely large enough (if necessary it will be enlarged)
  variant_subtree_nodes_expanded <- data.frame(matrix(data = 0, ncol = 5, nrow = max(2 * ceiling(max_cluster_size / detection_proba), ceiling(max_cluster_size / detection_proba) + 250)))
  names(variant_subtree_nodes_expanded) <- c("node_key", "mutation_occured", "variant_received", "variant_after_mutation", "detection")

  # create data frame to store information about the edges of the variant subtree we consider
  # we define it to be most likely large enough (if necessary it will be enlarged)
  variant_subtree_edges_expanded <- data.frame(matrix(data = 0, ncol = 5, nrow = max(2 * ceiling(max_cluster_size / detection_proba), ceiling(max_cluster_size / detection_proba) + 250)))
  names(variant_subtree_edges_expanded) <- c("from", "to", "variant_transmitted", "from_detected", "to_detected")


  mutation_proba <- 1 - exp(- yearly_mutation_rate / 365.25 * mean_generation_interval)
  detection_proba <- testing_proba * sequencing_proba

  variant_subtree_nodes_expanded[1,] <- c(1, stats::rbinom(n = 1, size = 1, prob = mutation_proba), 0, 1, stats::rbinom(n = 1, size = 1, prob = detection_proba))

  variant_subtree_free_leaves <- c(1)

  if (variant_subtree_nodes_expanded$mutation_occured[1] == 1) {

    variant <- 1

  } else {

    variant <- 0

  }

  # extend variant subtree ----

  # initialize the number of nodes, the number of free leaves and the number of edges of the variant subtree we consider
  number_nodes <- 1
  number_free_leaves <- 1
  number_edges <- 0

  # continue as long as there is a free leaf and the max_cluster_size limit of the variant subtree is not reached
  while ((number_free_leaves > 0) && (sum(variant_subtree_nodes_expanded$detection) <= max_cluster_size)) {

    # create new leaves (might be no new leaves)
    number_new_free_leaves <- stats::rnbinom(n = 1, size = k, mu = R)

    # add the new leaves to the set of nodes and to the set of free leaves and connect the new leaves to the graph
    if (number_new_free_leaves >= 1) {

      # write new leaves into a list
      new_leaves <- seq(from = max(variant_subtree_nodes_expanded$node_key) + 1, to = max(variant_subtree_nodes_expanded$node_key) + number_new_free_leaves, by = 1)

      # information about the newly generated nodes
      variant_subtree_nodes_expanded_temp <- data.frame(node_key = new_leaves,
                                                        mutation_occured = stats::rbinom(n = number_new_free_leaves, size = 1, prob = mutation_proba),
                                                        variant_received = rep(x = variant, times = number_new_free_leaves),
                                                        variant_after_mutation = new_leaves,
                                                        detection = stats::rbinom(n = number_new_free_leaves, size = 1, prob = detection_proba))

      # remove nodes at which a mutation occurred
      variant_subtree_nodes_expanded_temp <- variant_subtree_nodes_expanded_temp |> dplyr::filter(mutation_occured == 0)

      # update number of new free leaves (only those nodes at which no mutation occurred are considered)
      number_new_free_leaves <- nrow(variant_subtree_nodes_expanded_temp)

      if (number_new_free_leaves >= 1) {

        # update new_leaves
        new_leaves <- variant_subtree_nodes_expanded_temp |> dplyr::pull(node_key)

        # add variant_subtree_nodes_expanded_temp to variant_subtree_nodes_expanded
        variant_subtree_nodes_expanded[(number_nodes + 1):(number_nodes + number_new_free_leaves), ] <- variant_subtree_nodes_expanded_temp

        # add new free leaves to the set of nodes of the variant subtree and to the set of free leaves of the variant subtree
        variant_subtree_free_leaves[(number_free_leaves + 1):(number_free_leaves + number_new_free_leaves)] <- new_leaves

        # add edges connecting the new free leaves to the first free leaf of the tree
        variant_subtree_edges_expanded$from[(number_edges + 1):(number_edges + number_new_free_leaves)] <- rep(x = match(x = variant_subtree_free_leaves[1], table = variant_subtree_nodes_expanded$node_key), times = number_new_free_leaves)
        variant_subtree_edges_expanded$to[(number_edges + 1):(number_edges + number_new_free_leaves)] <- (number_edges + 1):(number_edges + number_new_free_leaves) + 1
        variant_subtree_edges_expanded$variant_transmitted[(number_edges + 1):(number_edges + number_new_free_leaves)] <- rep(x = variant, times = number_new_free_leaves)
        variant_subtree_edges_expanded$from_detected[(number_edges + 1):(number_edges + number_new_free_leaves)] <- unlist(lapply(X = variant_subtree_edges_expanded$from[(number_edges + 1):(number_edges + number_new_free_leaves)],
                                                                                                                                  FUN = function(x) if (variant_subtree_nodes_expanded$detection[x] == 1) {1} else {0}))
        variant_subtree_edges_expanded$to_detected[(number_edges + 1):(number_edges + number_new_free_leaves)] <- unlist(lapply(X = variant_subtree_edges_expanded$to[(number_edges + 1):(number_edges + number_new_free_leaves)],
                                                                                                                                FUN = function(x) if (variant_subtree_nodes_expanded$detection[x] == 1) {1} else {0}))

      }

      # update the number of nodes, the number of free leaves and the number of edges of the variant subtree
      number_nodes <- number_nodes + number_new_free_leaves
      number_free_leaves <- number_free_leaves + number_new_free_leaves
      number_edges <- number_edges + number_new_free_leaves

    }

    # delete the free leaf we are considering from the set of free leaves
    variant_subtree_free_leaves <- variant_subtree_free_leaves[-1]

    # update the number of free leaves of the tree
    number_free_leaves <- number_free_leaves - 1

  }

  # define output
  completed_variant_subtree <- list(nodes = variant_subtree_nodes_expanded |> dplyr::filter(node_key != 0),
                                    edges = variant_subtree_edges_expanded |> dplyr::filter(from != 0 & to != 0))

  return(completed_variant_subtree)

}


