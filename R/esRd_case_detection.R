
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom stats rbinom
esRd_case_detection_indep_bernoulli <- function(nodes, edges, detection_proba) {

  # make sure that the columns of nodes and edges have the right names
  names(nodes) <- c("node_key", "mutation_occured", "variant_received", "variant_after_mutation")
  names(edges) <- c("from", "to", "variant_transmitted")

  # prepare the data frames that will be returned
  nodes_a <- nodes |> dplyr::mutate(detection = rbinom(n = nrow(nodes), size = 1, prob = detection_proba))
  edges_a <- edges |> dplyr::mutate(from_detected = rep(x = 0, times = nrow(edges)),
                                    to_detected = rep(x = 0, times = nrow(edges)))

  # determine whether start and target of the edges have been detected
  if (nrow(edges_a) >= 1) {

    edges_a <- edges_a |> dplyr::mutate(from_detected = unlist(lapply(X = edges$from,
                                                                      FUN =  function(x) if (nodes_a |> dplyr::filter(node_key == x) |> dplyr::select(detection) == 1) {1} else {0})),
                                        to_detected = unlist(lapply(X = edges$to,
                                                                    FUN =  function(x) if (nodes_a |> dplyr::filter(node_key == x) |> dplyr::select(detection) == 1) {1} else {0})))

  }

  # create output
  tree_detection <- list(nodes = nodes_a, edges = edges_a)

  return(tree_detection)

}
