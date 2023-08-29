
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom stats rbinom
#' @importFrom stats rnbinom
esRd_complete_variant_subtree_old <- function(tree_nodes, tree_leaves, tree_edges, variant, limit_size, R, k, mutation_proba, detection_proba) {

  # create data frame to store information about the nodes of the variant subtree we consider
  # we define it to be large enough (CHECK: limit_size + 250, is there a better solution?)
  variant_subtree_nodes_expanded <- data.frame(matrix(data = 0, ncol = ncol(tree_nodes), nrow = limit_size + 250))
  names(variant_subtree_nodes_expanded) <- names(tree_nodes)

  # create data frame to store information about the edges of the variant subtree we consider
  # we define it to be large enough (CHECK: limit_size + 250, is there a better solution?)
  variant_subtree_edges_expanded <- data.frame(matrix(data = 0, ncol = ncol(tree_edges), nrow = limit_size + 250))
  names(variant_subtree_edges_expanded) <- names(tree_edges)

  # store information about the variant subtree we consider (part of variant subtree that has already been created)
  variant_subtree_nodes_0 <- tree_nodes |> dplyr::filter((mutation_occured == 1 & variant_after_mutation == variant) | (mutation_occured == 0 & variant_received == variant))
  variant_subtree_edges_0 <- tree_edges |> dplyr::filter((from %in% (variant_subtree_nodes_0 |> dplyr::pull(node_key))) & (to %in% (variant_subtree_nodes_0 |> dplyr::pull(node_key))))
  variant_subtree_edges_0 <- variant_subtree_edges_0 |> dplyr::mutate(from = unlist(lapply(X = from,
                                                                                           FUN = function(x) match(x = x, table = (variant_subtree_nodes_0 |> dplyr::pull(node_key))))),
                                                                      to = unlist(lapply(X = to,
                                                                                         FUN = function(x) match(x = x, table = (variant_subtree_nodes_0 |> dplyr::pull(node_key))))))

  variant_subtree_nodes_expanded[1:nrow(variant_subtree_nodes_0), ] <- variant_subtree_nodes_0

  if (nrow(variant_subtree_edges_0) >= 1) {

    variant_subtree_edges_expanded[1:nrow(variant_subtree_edges_0), ] <- variant_subtree_edges_0

  }

  # determine the nodes of the variant subtree we consider that are leaves
  variant_subtree_free_leaves <- intersect(variant_subtree_nodes_0 |> dplyr::pull(node_key), tree_leaves |> dplyr::pull(node_key))

  # extend variant subtree ----

  # initialize the number of nodes, the number of free leaves and the number of edges of the variant subtree we consider
  number_nodes <- nrow(variant_subtree_nodes_0)
  number_free_leaves <- length(variant_subtree_free_leaves)
  number_edges <- nrow(variant_subtree_edges_0)

  # continue as long as there is a free leaf and the limit_size limit of the variant subtree is not reached
  while ((number_free_leaves > 0) && (number_nodes <= limit_size)) {

    # create new leaves (might be no new leaves)
    number_new_free_leaves <- stats::rnbinom(n = 1, size = k, mu = R)

    # add the new leaves to the set of nodes and to the set of free leaves and connect the new leaves to the graph
    if (number_new_free_leaves >= 1) {

      # write new leaves into a list
      new_leaves <- seq(from = max(variant_subtree_nodes_expanded$node_key) + 1, to = max(variant_subtree_nodes_expanded$node_key) + number_new_free_leaves, by = 1)

      # information about the newly generated nodes
      variant_subtree_nodes_expanded_temp <- data.frame(node_key = new_leaves,
                                                        mutation_occured = rbinom(n = number_new_free_leaves, size = 1, prob = mutation_proba),
                                                        variant_received = rep(x = variant, times = number_new_free_leaves),
                                                        variant_after_mutation = new_leaves,
                                                        detection = rbinom(n = number_new_free_leaves, size = 1, prob = detection_proba))

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


