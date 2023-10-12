
#' @title Simulate transmission chain
#'
#' @description A transmission chain of maximal size `max_chain_size` is simulated.
#'
#' Viral transmission (`R` and `k`), viral mutation
#' (`yearly_mutation_rate` and `mean_generation_interval`) and case detection
#' (`testing_proba` and `sequencing_proba`) are included in the simulation.
#'
#' @param max_chain_size Maximal size of the simulated transmission chain.
#' @param R Effective reproduction number.
#' @param k Dispersion parameter.
#' @param yearly_mutation_rate Number of mutations of the viral genome per year.
#' @param mean_generation_interval Mean generation interval.
#' @param testing_proba Probability that a case is confirmed by a test.
#' @param sequencing_proba Probability that the genome of a confirmed case is sequenced.
#'
#' @return A list containing three data frames, `nodes` (information about cases), `edges` (information about transmission events)
#' and `free_leaves` (information about cases whose offspring has not been determined yet).
#'
#' The data frame `nodes` contains the columns `node_key` (id of node), `generation` (generation of nodes),
#' `mutation_occurred` (whether or not a mutation occurred), `variant_received` (variant received from parent node),
#' `variant_after_mutation` (variant present at node in case a mutation occurs),
#' `current_variant` (variant of virus present at node at time of detection) and `detection` (whether or not a case is detected).
#'
#' The data frame `edges` contains the columns `from` (index of infector node in `nodes`), `to` (index of infectee node in `nodes`),
#' `generation` (generation of infectee node), `from_variant` (variant of virus present at infector node at time of detection),
#' `to_variant` (variant of virus present at infectee node at time of detection), `variant_transmitted` (variant transmitted from infector to infectee),
#' `from_detected` (whether or not the infector node has been detected) and `to_detected` (whether or not the infectee node has been detected).
#'
#' The data frame `free_leaves` contains the columns `node_key` (id of node), `generation` (generation of nodes),
#' `mutation_occurred` (whether or not a mutation occurred), `variant_received` (variant received from parent node) and
#' `variant_after_mutation` (variant present after mutation at node in case a mutation occurs at that node).
#'
#' @export
#'
#' @examples estRodis_simulate_transmission_chain()

estRodis_simulate_transmission_chain <- function(max_chain_size = 100,
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

  # create data frames to store information about nodes (cases) and edges (transmissions) of the transmission chain
  # the data frames are defined to be sufficiently large, zeros are used as placeholders
  # we start with a single case
  nodes <- data.frame(node_key = c(1, rep(0, times = max_chain_size - 1)),
                      generation = c(1, rep(0, times = max_chain_size - 1)),
                      mutation_occurred = stats::rbinom(n = max_chain_size, size = 1, prob = mutation_proba),
                      variant_received = rep(0, times = max_chain_size),
                      variant_after_mutation = 1:max_chain_size,
                      current_variant = rep(0, times = max_chain_size),
                      detection = stats::rbinom(n = max_chain_size, size = 1, prob = detection_proba))

  edges <- data.frame(from = rep(x = 0, times = max_chain_size - 1),
                      to = rep(x = 0, times = max_chain_size - 1),
                      generation = rep(x = 0, times = max_chain_size - 1),
                      from_variant = rep(x = 0, times = max_chain_size - 1),
                      to_variant = rep(x = 0, times = max_chain_size - 1),
                      variant_transmitted = rep(x = 0, times = max_chain_size - 1),
                      from_detected = rep(x = 0, times = max_chain_size - 1),
                      to_detected = rep(x = 0, times = max_chain_size - 1))

  # in this list the nodes whose offspring have not yet been determined are stored
  free_leaves <- data.frame(node_key = c(1, rep(0, times = max_chain_size - 1)),
                            generation = c(1, rep(0, times = max_chain_size - 1)),
                            mutation_occurred = nodes$mutation_occurred,
                            variant_received = 0,
                            variant_after_mutation = 1:max_chain_size)

  # initialize the number of nodes, the number of free leaves and the number of edges of the transmission chain
  number_nodes <- 1
  number_free_leaves <- 1
  number_edges <- 0

  # continue as long as there is a free leaf and the size limit of the transmission chain (max_chain_size nodes) is not reached
  while ((number_free_leaves > 0) && (number_nodes <= max_chain_size)) {

    # create new leaves (might be no new leaves)
    number_new_free_leaves <- stats::rnbinom(n = 1, size = k, mu = R)

    if (number_new_free_leaves >= 1) {

      # write new leaves into a list, they are named iteratively
      new_leaves <- seq(from = number_nodes + 1, to = number_nodes + number_new_free_leaves, by = 1)

    }

    # add the new leaves to the set of nodes and to the set of free leaves and connect the new leaves to their parent node
    if (number_new_free_leaves != 0) {

      if (number_nodes + number_new_free_leaves <= max_chain_size) {

        # add new free leaves to the set of nodes of the transmission chain
        nodes$node_key[(number_nodes + 1):(number_nodes + number_new_free_leaves)] <- new_leaves
        nodes$generation[(number_nodes + 1):(number_nodes + number_new_free_leaves)] <- free_leaves$generation[1] + 1
        nodes$variant_received[(number_nodes + 1):(number_nodes + number_new_free_leaves)] <- (1 - free_leaves$mutation_occurred[1]) * free_leaves$variant_received[1] + free_leaves$mutation_occurred[1] * free_leaves$variant_after_mutation[1]

        # add new free leaves to the set of free leaves of the transmission chain
        free_leaves$node_key[(number_free_leaves + 1):(number_free_leaves + number_new_free_leaves)] <- new_leaves
        free_leaves$generation[(number_free_leaves + 1):(number_free_leaves + number_new_free_leaves)] <- free_leaves$generation[1] + 1
        free_leaves$variant_received[(number_free_leaves + 1):(number_free_leaves + number_new_free_leaves)] <- (1 - free_leaves$mutation_occurred[1]) * free_leaves$variant_received[1] + free_leaves$mutation_occurred[1] * free_leaves$variant_after_mutation[1]

        # add edges connecting the new free leaves to the first free leaf of the transmission chain
        edges$from[(number_edges + 1):(number_edges + number_new_free_leaves)] <- rep(x = free_leaves$node_key[1], times = number_new_free_leaves)
        edges$to[(number_edges + 1):(number_edges + number_new_free_leaves)] <- new_leaves
        edges$generation[(number_edges + 1):(number_edges + number_new_free_leaves)] <- free_leaves$generation[1] + 1

        # update the number of nodes, the number of free leaves and the number of edges of the transmission chain
        number_nodes <- number_nodes + number_new_free_leaves
        number_free_leaves <- number_free_leaves + number_new_free_leaves
        number_edges <- number_edges + number_new_free_leaves

        # delete the free leaf we are considering from the set of free leaves
        free_leaves <- free_leaves[-1,]

        # update the number of free leaves of the transmission chain
        number_free_leaves <- number_free_leaves - 1

      } else {

        # set number_nodes to a high enough value such that the condition for the while loop is no longer fulfilled
        number_nodes <- max_chain_size + 1

      }

    } else {

      # delete the free leaf we are considering from the set of free leaves
      free_leaves <- free_leaves[-1,]

      # update the number of free leaves of the tree
      number_free_leaves <- number_free_leaves - 1

    }

  }

  # remove rows containing zero from the data frames nodes, edges and free_leaves
  nodes <- nodes |> dplyr::filter(nodes$node_key != 0)

  edges <- edges |> dplyr::filter(edges$from != 0 & edges$to != 0)

  free_leaves <- free_leaves |> dplyr::filter(free_leaves$node_key != 0)


  # add column current_variant to nodes
  nodes <- nodes |> dplyr::mutate(current_variant = (1 - nodes$mutation_occurred) * nodes$variant_received + nodes$mutation_occurred * nodes$variant_after_mutation)

  if (nrow(edges) >= 1) {

    # add information about mutation process to edges
    edges <- edges |> dplyr::mutate(from_variant = nodes$current_variant[edges$from],
                                    to_variant = nodes$current_variant[edges$to],
                                    variant_transmitted = nodes$variant_received[edges$to])

    # add information about detection process to edges
    edges <- edges |> dplyr::mutate(from_detected = nodes$detection[edges$from],
                                    to_detected = nodes$detection[edges$to])

  }

  # create output
  structure_transmission_chain <- list(nodes = nodes,
                                       edges = edges,
                                       free_leaves = free_leaves)

  return(structure_transmission_chain)

}
