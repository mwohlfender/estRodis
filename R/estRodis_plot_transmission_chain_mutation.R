
# Warning messages:
# 1: Using size for a discrete variable is not advised.
# 2: Using the `size` aesthetic in this geom was deprecated in ggplot2 3.4.0.
# ℹ Please use `linewidth` in the `default_aes` field and elsewhere instead.
# This warning is displayed once every 8 hours.
# Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

#' @title Plot transmission chain with mutation
#'
#' @description A plot of a transmission chain including the mutation process is created.
#'
#' @param transmission_chain_nodes Data frame containing information about the nodes of the transmission chain
#' @param transmission_chain_edges Data frame containing information about the edges of the transmission chain
#' @param max_generation Maximum generation of nodes and edges that are included into the plot
#' @param style_plot Positioning of transmission chain on canvas
#' @param style_legend_clusters Level of detail of the legend of the identical sequence clusters
#'
#' @details The data frame `transmission_chain_nodes` has to contain at least the columns `node_key` (id of node),
#' `generation` (generation of nodes), `mutation_occurred` (whether or not a mutation occurred),
#' `variant_received` (variant received from parent node),
#' `variant_after_mutation` (variant present at node in case a mutation occurs) and
#' `current_variant` (variant of virus present at node at time of detection).
#'
#' The data frame `transmission_chain_edges` has to contain at least the columns
#' `from` (index of infector node in `transmission_chain_nodes`),
#' `to` (index of infectee node in `transmission_chain_nodes`),
#' `generation` (generation of infectee node),
#' `from_variant` (variant of virus present at infector node at time of detection),
#' `to_variant` (variant of virus present at infectee node at time of detection), and
#' `variant_transmitted` (variant transmitted from infector to infectee).
#'
#' Using the inputs `max_generation`, `style_plot` and `style_legend_clusters`,
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
#' Setting `style_legend_clusters = fixed`, the legend of the
#' identical sequence clusters includes all identical sequence clusters
#' of the transmission chain defined by nodes and edges regardless of their generation.
#'
#' Setting `style_legend_clusters = flexible`, the legend of the
#' identical sequence clusters only includes all identical sequence clusters
#' of the transmission chain defined by nodes and edges up to (and including)
#' generation `max_generation`.
#'
#' @return Plot of the transmission chain defined by nodes, edges and max_generation in the form of a cowplot plot_grid object.
#' @export

estRodis_plot_transmission_chain_mutation <- function(transmission_chain_nodes,
                                                      transmission_chain_edges,
                                                      max_generation = max(transmission_chain_nodes %>% dplyr::pull("generation")),
                                                      style_plot = "flexible",
                                                      style_legend_clusters = "flexible") {

  # #' @examples transmission_chain <- estRodis_simulate_transmission_chain()
  # #' transmission_chain_nodes <- transmission_chain$nodes
  # #' transmission_chain_edges <- transmission_chain$edges
  # #' plot_transmission_chain <- estRodis_plot_transmission_chain_mutation(
  # #' transmission_chain_nodes = transmission_chain_nodes,
  # #' transmission_chain_edges = transmission_chain_edges,
  # #' max_generation = min(max(transmission_chain_nodes$generation), 5))

  # determine nodes and edges that will be added to the plot
  if (style_plot == "fixed") {

    # the whole transmission chain (determined by `transmission_chain_nodes` and `transmission_chain_edges`) will be added to the plot,
    # but only nodes and edges up to generation `max_generation` will be visible
    plot_nodes <- transmission_chain_nodes
    plot_edges <- transmission_chain_edges

  } else {

    # only nodes and edges up to generation `max_generation` will be added to the plot
    plot_nodes <- transmission_chain_nodes %>% dplyr::filter(transmission_chain_nodes$generation <= max_generation)
    plot_edges <- transmission_chain_edges %>% dplyr::filter(transmission_chain_edges$generation <= max_generation)

  }

  # determine colors of the identical sequence clusters
  if (style_legend_clusters == "fixed") {

    # colors of all clusters of the whole transmission chain (determined by `transmission_chain_nodes` and `transmission_chain_edges`) will be added to the legend
    colors_clusters <- paletteer::paletteer_c("viridis::plasma", n = length(unique(transmission_chain_nodes %>% dplyr::pull("current_variant"))))

  } else {

    # only colors of the clusters visible on the plot (nodes and edges up to generation `max_generation`) will be added to the legend
    colors_clusters <- paletteer::paletteer_c("viridis::plasma", n = length(unique(transmission_chain_nodes %>% dplyr::filter(transmission_chain_nodes$generation <= max_generation) %>% dplyr::pull("current_variant"))))

  }

  # define colors for edges visible on the plot (edges up to generation max_generation)
  colors_transmission_chain_edges <- colors_clusters[unique(match(x = transmission_chain_edges %>% dplyr::filter(transmission_chain_edges$generation <= max_generation) %>% dplyr::pull("variant_transmitted"),
                                                                  table = unique(transmission_chain_nodes %>% dplyr::filter(transmission_chain_nodes$generation <= max_generation) %>% dplyr::pull("current_variant"))))]

  # define graph
  graph_transmission_chain <- tidygraph::tbl_graph(nodes = plot_nodes,
                                                   edges = plot_edges,
                                                   directed = TRUE)

  if (nrow(plot_edges) >= 1) {

    # add column "edge.id"
    plot_edges <- plot_edges %>% dplyr::mutate(edge.id = 1:nrow(plot_edges))

    # add x and y coordinate of nodes in column "from"
    plot_edges <- plot_edges %>% dplyr::mutate(node_key = plot_edges$from) %>% dplyr::left_join(create_layout(graph_transmission_chain, layout = "tree")  %>% dplyr::select(c("node_key", "x", "y")),
                                                                                              by = "node_key") %>% dplyr::rename("x_from" = "x", "y_from" = "y") %>% dplyr::select(-"node_key")

    # add x and y coordinate of nodes in column "to"
    plot_edges <- plot_edges %>% dplyr::mutate(node_key = plot_edges$to) %>% dplyr::left_join(create_layout(graph_transmission_chain, layout = "tree")  %>% dplyr::select(c("node_key", "x", "y")),
                                                                                            by = "node_key") %>% dplyr::rename("x_to" = "x", "y_to" = "y") %>% dplyr::select(-"node_key")

    # add column "in_same_cluster": 1 if both nodes of an edge belong to the same identical sequence cluster, 0 if not
    plot_edges <- plot_edges %>% dplyr::mutate(in_same_cluster = as.numeric(plot_edges$from_variant == plot_edges$to_variant))

  }

  # define breaks, limits, values and labels for the color and fill scales of the plot
  if (style_plot == "fixed" & style_legend_clusters == "fixed") {

    scale_color_fill_manual_breaks <- unique(c(transmission_chain_nodes$current_variant , -1))
    scale_color_fill_manual_limits <- factor(scale_color_fill_manual_breaks)
    scale_color_fill_manual_values <- c(colors_clusters, "white")
    scale_color_fill_manual_labels <- c(LETTERS[seq(from = 1, to = length(colors_clusters))], "")

  }

  if (style_plot == "fixed" & style_legend_clusters != "fixed") {

    scale_color_fill_manual_breaks <- unique(c(transmission_chain_nodes$current_variant * (transmission_chain_nodes$generation <=  max_generation) - (transmission_chain_nodes$generation > max_generation), -1))
    scale_color_fill_manual_limits <- factor(scale_color_fill_manual_breaks)
    scale_color_fill_manual_values <- c(colors_clusters, "white")
    scale_color_fill_manual_labels <- c(LETTERS[seq(from = 1, to = length(colors_clusters))], "")

  }

  if (style_plot != "fixed" & style_legend_clusters == "fixed") {

    scale_color_fill_manual_breaks <- unique(transmission_chain_nodes$current_variant)
    scale_color_fill_manual_limits <- factor(scale_color_fill_manual_breaks)
    scale_color_fill_manual_values <- colors_clusters
    scale_color_fill_manual_labels <- c(LETTERS[seq(from = 1, to = length(colors_clusters))])

  }

  if (style_plot != "fixed" & style_legend_clusters != "fixed") {

    scale_color_fill_manual_breaks <- unique(plot_nodes$current_variant)
    scale_color_fill_manual_limits <- factor(scale_color_fill_manual_breaks)
    scale_color_fill_manual_values <- colors_clusters
    scale_color_fill_manual_labels <- c(LETTERS[seq(from = 1, to = length(colors_clusters))])

  }

  # create plot
  plot_graph_mutation <- ggraph::ggraph(graph = graph_transmission_chain, layout = "tree") +
    {if (nrow(plot_edges) >= 1)
      ggraph::geom_edge_link(data = plot_edges,
                             mapping = ggplot2::aes(x = .data$x_from,
                                                    y = .data$y_from,
                                                    xend = .data$x_to,
                                                    yend = .data$y_to,
                                                    edge_color = factor(.data$variant_transmitted * (.data$generation <= max_generation) - (.data$generation > max_generation)),
                                                    edge_linetype = factor(.data$in_same_cluster),
                                                    edge_width = factor(.data$in_same_cluster)))} +
    {if (nrow(plot_edges) >= 1)
      ggraph::scale_edge_color_manual(name = "Identical sequence cluster:",
                                      breaks = unique(c(plot_edges$variant_transmitted * (plot_edges$generation <= max_generation) - (plot_edges$generation > max_generation), -1)),
                                      limits = factor(unique(c(plot_edges$variant_transmitted * (plot_edges$generation <= max_generation) - (plot_edges$generation > max_generation), -1))),
                                      values = c(colors_transmission_chain_edges, "white"),
                                      guide = "none")} +
    {if (nrow(plot_edges) >= 1)
      ggraph::scale_edge_width_manual(name = "Transmission",
                                      breaks = c(1, 0),
                                      limits = factor(c(1, 0)),
                                      values = c(1, 0.5),
                                      guide = "none")} +
    {if (nrow(plot_edges) >= 1)
      ggraph::scale_edge_linetype_manual(name = "Transmission:",
                                         breaks = c(1, 0),
                                         limits = factor(c(1, 0)),
                                         values = c("solid", "dashed"),
                                         labels = c("Within same cluster", "Not within same cluster"))} +
    ggraph::geom_node_point(data = create_layout(graph_transmission_chain, layout = "tree"),
                            mapping = ggplot2::aes(x = .data$x,
                                                   y = .data$y,
                                                   color = factor(.data$current_variant * (.data$generation <=  max_generation) - (.data$generation > max_generation)),
                                                   fill = factor(.data$current_variant * (.data$generation <= max_generation) - (.data$generation > max_generation)),
                                                   shape = factor(.data$mutation_occurred),
                                                   size = factor(.data$mutation_occurred))) +
    ggplot2::scale_color_manual(name = "Identical sequence cluster:",
                                breaks = scale_color_fill_manual_breaks,
                                limits = scale_color_fill_manual_limits,
                                values = scale_color_fill_manual_values,
                                labels = scale_color_fill_manual_labels) +
    ggplot2::scale_fill_manual(name = "Cluster",
                               breaks = scale_color_fill_manual_breaks,
                               limits = scale_color_fill_manual_limits,
                               values = scale_color_fill_manual_values,
                               labels = scale_color_fill_manual_labels,
                               guide = "none") +
    ggplot2::scale_shape_manual(name = "Mutation:",
                                breaks = c(1, 0),
                                limits = factor(c(1, 0)),
                                values = c(22, 21),
                                labels = c("Mutation occurred", "No mutation occurred")) +
    ggplot2::scale_size_manual(name = "Mutation:",
                               breaks = c(1, 0),
                               limits = factor(c(1, 0)),
                               values = c(3.5, 2.5),
                               guide = "none") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = "black"),
                   legend.position = "right",
                   legend.box = "vertical",
                   legend.margin = ggplot2::margin(t = 0.01, r = 0.05, b = 0.01, l = 0.05, unit = "cm"),
                   legend.key = ggplot2::element_rect(fill = "white")) +
    ggplot2::guides(color = "none",
                    edge_linetype = ggplot2::guide_legend(order = 1),
                    shape = ggplot2::guide_legend(order = 2, override.aes=list(fill="black")))


  plot_graph_mutation_1 <- plot_graph_mutation +
    ggplot2::guides(color = ggplot2::guide_legend(order = 1, nrow = 1),
                    edge_linetype = "none",
                    shape = "none") +
    ggplot2::theme(legend.box = "vertical",
                   legend.margin = ggplot2::margin(t = 0, r = 0.1, b = 0, l = 0.1, unit = "cm"),
                   legend.key = ggplot2::element_rect(fill = "white")) +
    ggplot2::theme(legend.spacing.x = ggplot2::unit(0.05, "cm"))


  plot_graph_mutation_guide_color <- cowplot::get_legend(plot_graph_mutation_1)


  output <- cowplot::plot_grid(plot_graph_mutation,
                               plot_graph_mutation_guide_color,
                               nrow = 2,
                               rel_heights = c(5, 1)) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"))


  return(output)

}
