#' Weighted bond partners around central node
#'
#' This function calculates the weighted bond partners around a central node in
#' a given \code{\link[tidygraph]{tbl_graph}} object. The weighted bond partners
#' are determined based on a specified sphere around the central node. The
#' resulting data frame includes the symbols (atom types) and the weighted
#' distances of the bond partners.
#' 
#' @param tbl_graph The input tbl_graph object representing the molecular graph.
#' @param center An integer specifying the central node.
#' @param sphere An integer indicating how many nodes around the center should
#'   be included.
#' 
#' @return Returns a \code{data.frame} containing the symbols (atom types) and the
#' edge-weighted distances of the bond partners.
#' 
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
#'   gr <- sdf2tbl_graph(molec)
#'   wbp <- weight_bond_partners(gr, center = 3, sphere = 1)
#' }
#' 
#' @export
#' 
#' @importFrom tidygraph activate
#' @importFrom igraph distances
#' @importFrom magrittr "%>%"
#' @importFrom dplyr rename mutate select arrange
weight_bond_partners <- function(tbl_graph, center = 1, sphere = 1) {
  d <- NULL # To avoid NOTE "no visible binding for global variable"
  .G <- getFromNamespace(".G", "tidygraph")
  res <- tbl_graph %>%
    activate("edges") %>%
    rename(weight = "bondorder") %>%
    activate("nodes") %>%
    mutate(wd = distances(.G(), to = center)) %>%
    mutate(d = distances(.G(), to = center, weights = NA)) %>%
    as.data.frame() %>%
    subset(d == sphere) %>%
    select("sym", "wd") %>%
    arrange("sym", "wd")
  
  return(res)
}

#' Split molecular graph by removing a specified atom type
#'
#' This function splits the given \code{\link[tidygraph]{tbl_graph}} object by
#' removing all by-type specified atoms from the graph. The resulting graph(s)
#' consist(s) of (multiple) connected component(s), each representing a separate
#' molecular graph after the split.
#' 
#' @param tbl_graph The input tbl_graph object representing the molecular graph.
#' @param splitpos A character specifying the atom type to be removed by the
#'   splitting process.
#'
#' @return Returns a list of tbl_graph objects, each representing a separate
#'   molecular graph after the split.
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
#'   gr <- sdf2tbl_graph(molec)
#'   subs <- split_graph_at(gr, splitpos = "P")
#' }
#' 
#' @export
#' 
#' @importFrom tidygraph morph crystallize to_components
#' @importFrom dplyr pull filter
#' @importFrom magrittr "%>%"
split_graph_at <- function(tbl_graph, splitpos = "P") {
  sym <- NULL # To avoid NOTE "no visible binding for global variable"
  
  resu <- tbl_graph %>%
    filter(sym != splitpos) %>%
    morph(to_components, type = "strong") %>%
    crystallize() %>%
    pull("graph")
  
  return(resu)
}

#' Subset molecular graph by focusing on the neighborhood around one atom
#'
#' This function subsets the given \code{\link[tidygraph]{tbl_graph}} object by
#' focusing on the neighborhood around a specific atom. The resulting graph
#' includes only the atoms within a certain distance from the specified atom.
#' 
#' @param tbl_graph The input tbl_graph object representing the molecular graph.
#' @param focuspos A character specifying the atom type to focus on.
#' @param loc_neigh A positive integer indicating how many bonds away an atom
#'   can be from the central atom to be considered.
#'
#' @return Returns a tbl_graph object representing the subset of the molecular
#'   graph based on the neighborhood around the specified atom.
#'   
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
#'   gr <- sdf2tbl_graph(molec)
#'   subs <- subset_graph_at(gr, focuspos = "P", loc_neigh = 1)
#' }
#' 
#' @export
#' 
#' @importFrom tidygraph morph crystallize to_local_neighborhood
#' @importFrom dplyr pull filter
#' @importFrom magrittr "%>%" extract2
subset_graph_at <- function(tbl_graph, focuspos = "P", loc_neigh = 3) {
  sym <- NULL # To avoid NOTE "no visible binding for global variable"
  
  s_pos <- tbl_graph %>%
    filter(sym == focuspos) %>%
    select("id") %>%
    data.frame() %>% as.numeric()
  
  resu <- tbl_graph %>%
    morph(to_local_neighborhood, s_pos, order = loc_neigh, mode = "all") %>%
    crystallize() %>%
    pull("graph") %>%
    extract2(1)
  
  return(resu)
}