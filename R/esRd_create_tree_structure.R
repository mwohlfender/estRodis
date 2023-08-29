
#' @importFrom dplyr filter
#' @importFrom stats rnbinom
esRd_create_tree_structure <- function(R, k, max_tree_size) {

  # create initial graph (just one node) (data frames are defined to be sufficiently large, zeros are used as placeholders)
  nodes_tree <- c(1, rep(0, times = max_tree_size - 1))
  free_leaves_tree <- c(1, rep(0, times = max_tree_size - 1))
  edges_tree <- data.frame(from = rep(x = 0, times = max_tree_size - 1), to = rep(x = 0, times = max_tree_size - 1))

  # initialize the number of nodes, the number of free leaves and the number of edges of the tree
  number_nodes <- 1
  number_free_leaves <- 1
  number_edges <- 0

  # continue as long as there is a free leaf and the size limit of the tree (max_tree_size nodes) is not reached
  while ((number_free_leaves > 0) && (number_nodes <= max_tree_size)) {

    # create new leaves (might be no new leaves)
    number_new_free_leaves <- stats::rnbinom(n = 1, size = k, mu = R)

    if (number_new_free_leaves >= 1) {

      # write new leaves into a list
      new_leaves <- seq(from = number_nodes + 1, to = number_nodes + number_new_free_leaves, by = 1)

    }

    # add the new leaves to the set of nodes and to the set of free leaves and connect the new leaves to the graph
    if (number_new_free_leaves != 0) {

      if (number_nodes + number_new_free_leaves <= max_tree_size) {

        # add new free leaves to the set of nodes of the tree and to the set of free leaves of the tree
        nodes_tree[(number_nodes + 1):(number_nodes + number_new_free_leaves)] <- new_leaves
        free_leaves_tree[(number_free_leaves + 1):(number_free_leaves + number_new_free_leaves)] <- new_leaves

        # add edges connecting the new free leaves to the first free leaf of the tree
        edges_tree$from[(number_edges + 1):(number_edges + number_new_free_leaves)] <- rep(x = free_leaves_tree[1], times = number_new_free_leaves)
        edges_tree$to[(number_edges + 1):(number_edges + number_new_free_leaves)] <- new_leaves

        # update the number of nodes, the number of free leaves and the number of edges of the tree
        number_nodes <- number_nodes + number_new_free_leaves
        number_free_leaves <- number_free_leaves + number_new_free_leaves
        number_edges <- number_edges + number_new_free_leaves

        # delete the free leaf we are considering from the set of free leaves
        free_leaves_tree <- free_leaves_tree[-1]

        # update the number of free leaves of the tree
        number_free_leaves <- number_free_leaves - 1

      } else {

        # set number_nodes to a high enough value such that the condition for the while loop is no longer fulfilled
        number_nodes <- max_tree_size + 1

      }

    } else {

      # delete the free leaf we are considering from the set of free leaves
      free_leaves_tree <- free_leaves_tree[-1]

      # update the number of free leaves of the tree
      number_free_leaves <- number_free_leaves - 1

    }

  }

  # create output: remove rows containing zero from the data frames nodes_tree, free_leaves_tree and edges_tree
  tree_structure <- list(nodes = data.frame(node_key = nodes_tree[nodes_tree != 0]),
                         free_leaves = data.frame(node_key = free_leaves_tree[free_leaves_tree != 0]),
                         edges = edges_tree |> dplyr::filter(from != 0 & to != 0))

  return(tree_structure)

}
