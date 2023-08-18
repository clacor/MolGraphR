#' Plot molecular tbl_graph using ggraph
#'
#' This function plots a molecular \code{\link[tidygraph]{tbl_graph}} using the
#' \code{\link[ggraph]{ggraph}} package. The function generates a visualization
#' of the molecular graph, where nodes represent atoms and edges represent
#' bonds. The appearance of the plot can be customized using the provided color
#' scheme and van der Waals radii.
#' 
#' @param tbl_graph The input tbl_graph object representing the molecular graph.
#' @param my_colors A vector of colors to assign to different atom types. The
#'   default color scheme includes predefined colors for common atom types H, C,
#'   P, O, N, S, F, Cl, and I.
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
#'   gr <- sdf2tbl_graph(molec)
#'   mol_plot(gr)
#' }
#' 
#' @export
#' 
#' @importFrom ggraph ggraph geom_edge_link scale_edge_width geom_node_point theme_graph
#' @importFrom ggplot2 scale_color_manual labs aes
#' @importFrom PeriodicTable rvdw
mol_plot <- function(tbl_graph, my_colors = c(H = "grey80", C = "grey30",
                                              P = "blue", O = "red", N = "green", S = "orange", F = "cyan",
                                              Cl = "chartreuse", I = "darkmagenta")) {
  sym <- bondorder <- NULL # To avoid NOTE "no visible binding for global variable"
  radii <- PeriodicTable::rvdw(tbl_graph %>% activate("nodes") %>% data.frame() %>% pull("sym"))
  ggraph(tbl_graph, layout = "stress") +
    geom_edge_link(aes(width = bondorder), alpha = 0.8) +
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_point(aes(col = sym, size = radii)) +
    scale_color_manual(values = my_colors) +
    theme_graph() +
    labs(col = "Atom Type", size = "Van der Waals Radius", width = "Bond Order")
}

#' Plot molecular tbl_graph using ggraph
#'
#' This function plots a molecular \code{\link[tidygraph]{tbl_graph}} using the
#' \code{\link[ggraph]{ggraph}} package, focusing on a central node. The
#' resulting visualization shows the molecular graph with the central node at
#' the center, and the distances from the central node are reflected in the node
#' sizes. The appearance of the plot can be customized using the provided color
#' scheme.
#' 
#' @param tbl_graph The input tbl_graph object representing the molecular graph.
#' @param center An integer specifying the central node to focus on.
#' @param my_colors A vector of colors to assign to different atom types. The
#'   default color scheme includes predefined colors for common atom types H, C,
#'   P, O, N, S, F, Cl, and I.
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
#'   gr <- sdf2tbl_graph(molec)
#'   mol_plot_center(gr, center = 3)
#' }
#'
#' @export
#' 
#' @importFrom utils getFromNamespace
#' @importFrom ggraph ggraph geom_edge_link scale_edge_width geom_node_point theme_graph
#' @importFrom ggplot2 scale_color_manual labs coord_fixed aes
#' @importFrom ggforce geom_circle
#' @importFrom igraph distances
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate
mol_plot_center <- function(tbl_graph, center = 1, my_colors = c(H = "grey80",
                                                                 C = "grey30", P = "blue", O = "red", N = "green", S = "orange", F = "cyan",
                                                                 Cl = "chartreuse", I = "darkmagenta")) {
  d <- r <- sym <- bondorder <- NULL # To avoid NOTE "no visible binding for global variable"
  .G <- getFromNamespace(".G", "tidygraph")
  tbl_graph %>%
    mutate(d = distances(.G(), to = center)) %>%
    ggraph(layout = "focus", focus = center) +
    geom_edge_link(aes(width = bondorder), alpha = 0.8) +
    scale_edge_width(range = c(0.2, 2), breaks = c(1, 2, 3)) +
    geom_circle(aes(x0 = 0, y0 = 0, r = r), data.frame(r = 1:3), colour = 'grey') + 
    geom_node_point(aes(col = sym, size = 3*1/(d+1))) +
    coord_fixed() + 
    scale_color_manual(values = my_colors) +
    theme_graph() +
    labs(col = "Atom Type", size = "Distance from Center", edge_width = "Bond Order")
}