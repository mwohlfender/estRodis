

#' Plot transmission chain
#'
#' @description A plot of a transmission chain is created.
#'
#' @param nodes Data frame containing information about the nodes of the transmission chain
#' @param edges Data frame containing information about the edges of the transmission chain
#' @param max_generation Maximum generation of nodes and edges that should be included into the plot
#' @param colour_nodes_edges Colour for nodes and edges
#'
#' @details The data frame `nodes` has to contain at least the columns `node_key` (id of node) and `generation` (generation of nodes).
#'
#' The data frame `edges` has to contain at least the columns `from` (index of infector node in `nodes`), `to` (index of infectee node in `nodes`),
#' `generation` (generation of infectee node),
#'
#' @return Plot of the transmission chain defined by nodes, edges and max_generation in the form of a ggplot2 object.
#' @export
#'
#' @examples transmission_chain <- estRodis_simulate_transmission_chain()
#' nodes <- transmission_chain$nodes |> dplyr::select(c("node_key", "generation"))
#' edges <- transmission_chain$edges |> dplyr::select(c("from", "to", "generation"))
#' plot_transmission_chain <- estRodis_plot_transmission_chain(nodes = nodes, edges = edges, max_generation = min(max(nodes$generation), 5))
estRodis_plot_transmission_chain <- function(nodes, edges, max_generation, colour_nodes_edges = "navyblue") {

  # define graph
  graph_transmission_tree <- tidygraph::tbl_graph(nodes = nodes,
                                                  edges = edges,
                                                  directed = TRUE)

  if (nrow(edges) >= 1) {

    # add x and y coordinate of nodes in column "from"
    edges <- edges |> dplyr::mutate(node_key = from) |> dplyr::left_join(ggraph::create_layout(graph_transmission_tree, layout = "tree")  |> dplyr::select(c("node_key", "x", "y")),
                                                                         by = "node_key") |> dplyr::rename(x_from = x, y_from = y) |> dplyr::select(-"node_key")

    # add x and y coordinate of nodes in column "to"
    edges <- edges |> dplyr::mutate(node_key = to) |> dplyr::left_join(ggraph::create_layout(graph_transmission_tree, layout = "tree")  |> dplyr::select(c("node_key", "x", "y")),
                                                                       by = "node_key") |> dplyr::rename(x_to = x, y_to = y) |> dplyr::select(-"node_key")

    # add column "edge.id"
    edges <- edges |> dplyr::mutate(edge.id = 1:nrow(edges))

  }

  # create plot
  plot_graph_mutation_detection <- ggraph::ggraph(graph = graph_transmission_tree, layout = "tree") +
    {if (nrow(edges) >= 1) ggraph::geom_edge_link(data = edges,
                                                  mapping = ggplot2::aes(x = x_from,
                                                                         y = y_from,
                                                                         xend = x_to,
                                                                         yend = y_to,
                                                                         edge_colour = factor((generation <= max_generation) - (generation > max_generation))),
                                                  edge_linetype = "solid",
                                                  edge_width = 1) } +
    {if (nrow(edges) >= 1) ggraph::scale_edge_color_manual(name = "Identical sequence cluster:",
                                                           breaks = c(1, -1),
                                                           limits = factor(c(1, -1)),
                                                           values = c(colour_nodes_edges, "white"),
                                                           guide = "none") } +
    ggraph::geom_node_point(data = ggraph::create_layout(graph_transmission_tree, layout = "tree"),
                            mapping = ggplot2::aes(x = x,
                                                   y = y,
                                                   color = factor((generation <= max_generation) - (generation > max_generation)),
                                                   fill = factor((generation <= max_generation) - (generation > max_generation))),
                            shape = 21,
                            size = 2.5) +
    ggplot2::scale_color_manual(name = "Identical sequence cluster:",
                                breaks = c(1, -1),
                                limits = factor(c(1, -1)),
                                values = c(colour_nodes_edges, "white"),
                                guide = "none") +
    ggplot2::scale_fill_manual(name = "Cluster",
                               breaks = c(1, -1),
                               values = c(colour_nodes_edges, "white"),
                               guide = "none") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", colour = "black"))

  return(plot_graph_mutation_detection)

}

# transmission_chain <- estRodis_simulate_transmission_chain()
# nodes <- transmission_chain$nodes |> dplyr::select(c("node_key", "generation"))
# edges <- transmission_chain$edges |> dplyr::select(c("from", "to", "generation"))
# plot_transmission_chain <- estRodis_plot_transmission_chain(nodes = nodes, edges = edges, max_generation = min(max(nodes$generation), 5))
# plot_transmission_chain
