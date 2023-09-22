#' The 'estRodis' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE
#'
#' @docType package
#' @name estRodis-package
#' @aliases estRodis
#' @useDynLib estRodis, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import RcppParallel
#' @import rstantools
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr select
#' @importFrom dplyr rename
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 theme
#' @importFrom ggraph create_layout
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_link
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph scale_edge_color_manual
#' @importFrom parallelly availableCores
#' @importFrom rlang .data
#' @importFrom rstan sampling
#' @importFrom stats rbinom
#' @importFrom stats rnbinom
#' @importFrom tidygraph tbl_graph
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.26.16. https://mc-stan.org
#'
NULL
