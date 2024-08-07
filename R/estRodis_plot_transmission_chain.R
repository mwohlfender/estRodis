
#' @title Plot transmission chain
#'
#' @description A plot of a transmission chain is created.
#'
#' @param transmission_chain_nodes Data frame containing information about the nodes of the transmission chain
#' @param transmission_chain_edges Data frame containing information about the edges of the transmission chain
#' @param max_generation Maximum generation of nodes and edges that are included into the plot
#' @param colour_nodes_edges Colour for nodes and edges
#' @param style_plot Positioning of transmission chain on canvas
#' @param width_edges Width of edges of graph
#'
#' @details The data frame `transmission_chain_nodes` has to contain at least the columns `node_key` (id of node)
#' and `generation` (generation of nodes).
#'
#' The data frame `transmission_chain_edges` has to contain at least the columns `from`
#' (index of infector node in `transmission_chain_nodes`),
#' `to` (index of infectee node in `transmission_chain_nodes`)
#' and `generation` (generation of infectee node).
#'
#' Using the inputs `max_generation` and `style_plot`,
#' the graphical presentation of the transmission chain can be influenced.
#'
#' `max_generation` determines the last generation of nodes and edges
#' that are included into the plot.
#'
#' Setting `style_plot = fixed`, nodes and edges up to (and including)
#' generation `max_generation` are drawn,
#' while nodes and edges of later generations are made invisible.
#'
#' Setting `style_plot = flexible`, nodes and edges up to (and including)
#' generation `max_generation` are drawn,
#' while nodes and edges of later generations are cut away.
#'
#' @return Plot of the transmission chain defined by nodes, edges and max_generation in the form of a ggplot2 object.
#' @export

estRodis_plot_transmission_chain <- function(transmission_chain_nodes,
                                             transmission_chain_edges,
                                             max_generation = max(transmission_chain_nodes |> dplyr::pull("generation")),
                                             colour_nodes_edges = "navyblue",
                                             style_plot = "flexible",
                                             width_edges = 1) {

  # #' @examples transmission_chain <- estRodis_simulate_transmission_chain()
  # #' transmission_chain_nodes <- transmission_chain$nodes |> dplyr::select(c("node_key", "generation"))
  # #' transmission_chain_edges <- transmission_chain$edges |> dplyr::select(c("from", "to", "generation"))
  # #' plot_transmission_chain <- estRodis_plot_transmission_chain(
  # #' transmission_chain_nodes = transmission_chain_nodes,
  # #' transmission_chain_edges = transmission_chain_edges,
  # #' max_generation = min(max(nodes$generation), 5))

  # determine nodes and edges that will be added to the plot
  if (style_plot == "fixed") {

    # the whole transmission chain (determined by `transmission_chain_nodes` and `transmission_chain_edges`) will be added to the plot,
    # but only nodes and edges up to generation `max_generation` will be visible
    plot_nodes <- transmission_chain_nodes
    plot_edges <- transmission_chain_edges

  } else {

    # only nodes and edges up to generation `max_generation` will be added to the plot
    plot_nodes <- transmission_chain_nodes |> dplyr::filter(transmission_chain_nodes$generation <= max_generation)
    plot_edges <- transmission_chain_edges |> dplyr::filter(transmission_chain_edges$generation <= max_generation)

  }

  # define graph
  graph_transmission_chain <- tidygraph::tbl_graph(nodes = plot_nodes,
                                                   edges = plot_edges,
                                                   directed = TRUE)

  if (nrow(plot_edges) >= 1) {

    # add column "edge.id"
    plot_edges <- plot_edges |> dplyr::mutate(edge.id = 1:nrow(plot_edges))

    # add x and y coordinate of nodes in column "from"
    plot_edges <- plot_edges |> dplyr::mutate(node_key = plot_edges$from) |> dplyr::left_join(create_layout(graph_transmission_chain, layout = "tree")  |> dplyr::select(c("node_key", "x", "y")),
                                                                                              by = "node_key") |> dplyr::rename("x_from" = "x", "y_from" = "y") |> dplyr::select(-"node_key")

    # add x and y coordinate of nodes in column "to"
    plot_edges <- plot_edges |> dplyr::mutate(node_key = plot_edges$to) |> dplyr::left_join(create_layout(graph_transmission_chain, layout = "tree")  |> dplyr::select(c("node_key", "x", "y")),
                                                                                            by = "node_key") |> dplyr::rename("x_to" = "x", "y_to" = "y") |> dplyr::select(-"node_key")

  }


  # create plot
  plot_graph <- ggraph::ggraph(graph = graph_transmission_chain, layout = "tree") +
    {if (nrow(plot_edges) >= 1) ggraph::geom_edge_link(data = plot_edges,
                                                       mapping = ggplot2::aes(x = .data$x_from,
                                                                              y = .data$y_from,
                                                                              xend = .data$x_to,
                                                                              yend = .data$y_to,
                                                                              edge_colour = factor((.data$generation <= max_generation) - (.data$generation > max_generation))),
                                                       edge_linetype = "solid",
                                                       edge_width = width_edges)} +
    {if (nrow(plot_edges) >= 1) ggraph::scale_edge_color_manual(name = "Identical sequence cluster:",
                                                                breaks = c(1, -1),
                                                                limits = factor(c(1, -1)),
                                                                values = c(colour_nodes_edges, "white"),
                                                                guide = "none")} +
    ggraph::geom_node_point(data = ggraph::create_layout(graph_transmission_chain, layout = "tree"),
                            mapping = ggplot2::aes(x = .data$x,
                                                   y = .data$y,
                                                   color = factor((.data$generation <= max_generation) - (.data$generation > max_generation)),
                                                   fill = factor((.data$generation <= max_generation) - (.data$generation > max_generation))),
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



  return(plot_graph)

}


