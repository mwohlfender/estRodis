
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom stats rbinom
esRd_apply_mutations <- function(nodes, edges, mutation_proba) {

  # make sure that the columns of nodes and edges have the right names
  names(nodes) <- c("node_key")
  names(edges) <- c("from", "to")

  # determine the number of nodes
  number_nodes <- nrow(nodes)

  # determine values or add placeholders for the following properties of the nodes:
  # (a) whether a mutation occurs at the node,
  # (b) the variant the node received from its direct ancestor and
  # (c) the variant the node transmits to its direct offspring if a mutation occurs at the node
  nodes_a <- nodes |> dplyr::mutate(mutation_occured = rbinom(n = number_nodes, size = 1, prob = mutation_proba),
                                    variant_received = rep(x = 0, times = number_nodes),
                                    variant_after_mutation = nodes$node_key)

  # add placeholder for the variant transmitted via an edge
  edges_a <- edges |> dplyr::mutate(variant_transmitted = rep(x = 0, times = nrow(edges)))

  # if there is only one node, then there is not much to do
  if (number_nodes >= 2) {

    for (ii in 1:number_nodes) {

      # get edges connecting nodes_a$node_key[ii] to its direct offspring
      edges_to_direct_offspring <- edges |> dplyr::filter(from == ii)

      # determine direct offspring of nodes_a$node_key[ii]
      descendants_one <- unique(union(nodes$node_key[edges_to_direct_offspring$from], nodes$node_key[edges_to_direct_offspring$to]))
      descendants_one <- descendants_one[descendants_one != ii]

      # if nodes_a$node_key[ii] has no direct offspring, then there is nothing further to do
      if (length(descendants_one) >= 1) {

        for (jj in 1:length(descendants_one)) {

          # depending on whether a mutation takes place at a node or not, a different variant is transmitted to its direct offspring
          if (nodes_a$mutation_occured[ii]) {

            nodes_a$variant_received[descendants_one[jj]] <- nodes_a$variant_after_mutation[ii]
            edges_a$variant_transmitted[descendants_one[jj] - 1] <- nodes_a$variant_after_mutation[ii]

          } else {

            nodes_a$variant_received[descendants_one[jj]] <- nodes_a$variant_received[ii]
            edges_a$variant_transmitted[descendants_one[jj] - 1] <- nodes_a$variant_received[ii]

          }

        }

      }

    }

  }

  # create output
  tree_mutation <- list(nodes = nodes_a, edges = edges_a)

  return(tree_mutation)

}
