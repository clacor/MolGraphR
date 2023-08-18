#' Transform SDF format into a tbl_graph
#' 
#' Transforms \code{ChemmineR}'s \code{\link[ChemmineR]{SDFset}} format into a
#' \code{\link[tidygraph]{tbl_graph}} from \code{tidygraph}.
#'
#' This function takes a \code{ChemmineR} \code{\link[ChemmineR]{SDFset}}
#' (Structure-Data File) format and converts it into a
#' \code{\link[tidygraph]{tbl_graph}} object from the \code{tidygraph} package.
#' The resulting \code{\link[tidygraph]{tbl_graph}} represents a molecular graph
#' with nodes as atoms and edges as bonds.
#' 
#' @param sdf The input SDF object in ChemmineR format.
#' 
#' @return Returns a tbl_graph object representing the molecular graph.
#' 
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
#'   gr <- sdf2tbl_graph(molec)
#' }
#' 
#' @export
#' 
#' @importFrom stringr str_split
#' @importFrom tidygraph tbl_graph
#' @importFrom ChemmineR atomblock bondblock
sdf2tbl_graph <- function(sdf) {
  temp <- do.call(rbind, str_split(rownames(atomblock(sdf)[[1]]), "_"))
  nodes <- data.frame(id = as.numeric(temp[, 2]), sym = temp[, 1])
  
  edges <- data.frame(bondblock(sdf)[[1]][, 1:3])
  colnames(edges) <- c("A1", "A2", "bondorder")
  
  molgraph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  
  return(molgraph)
}

#' Transform a tbl_graph into SDF format
#' 
#' Transform a \code{\link[tidygraph]{tbl_graph}} from \code{tidygraph} to
#' \code{ChemmineR} \code{\link[ChemmineR]{SDFset}} format
#'
#' This function takes a \code{\link[tidygraph]{tbl_graph}} object from the
#' \code{tidygraph} package and converts it into \code{ChemmineR}'s
#' \code{\link[ChemmineR]{SDFset}} (Structure-Data File) format. The resulting
#' SDF object can be used for further analysis or storage.
#' 
#' @param tbl_graph The input tbl_graph object representing the molecular graph.
#' 
#' @return Returns an SDFset object in ChemmineR format.
#' 
#' @examples 
#' \dontrun{
#'   ## NOT RUN
#'   molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
#'   gr <- sdf2tbl_graph(molec)
#'   molec2 <- tbl_graph2sdf(gr)
#' }
#' 
#' @export
#' 
#' @importFrom methods new
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggraph ggraph geom_node_point
#' @importFrom tidygraph activate
#' @importFrom magrittr "%>%"
tbl_graph2sdf <- function(tbl_graph) {
  ## Geometry is not optimized
  xycoord <- ggplot_build(ggraph(tbl_graph, "stress") + geom_node_point())[[1]][[1]][, 1:2]

  AB <- tbl_graph %>% activate("nodes") %>% data.frame
  BB <- tbl_graph %>% activate("edges") %>% data.frame

  HL <- c(Molecule_Name = "", Source = "tbl_graph2sdf-from-R", Comment = "",
    Counts_Line = paste("", nrow(AB), nrow(BB), "0", "0", "0", "0", "0", "0",
                        "0", "0999", "V2000", sep = " "))

  ABnew <- as.matrix(cbind(xycoord, matrix(0, nrow(AB), 14)))
  colnames(ABnew) <- paste0("C", 1:16)
  rownames(ABnew) <- paste0(AB$sym, "_", 1:nrow(AB))

  BBnew <- as.matrix(cbind(BB[, 1:3], matrix(0, nrow(BB), 4)))
  colnames(BBnew) <- paste0("C", 1:7)
  rownames(BBnew) <- if(nrow(BBnew) == 0){NULL} else{1:nrow(BBnew)}

  sdf <- list(new("SDF", header = HL, atomblock = ABnew, bondblock = BBnew))

  sdfset <- new("SDFset", SDF = sdf)

  return(sdfset)
}

#' Get molecular formula from molecular tbl_graph
#' 
#' This function calculates the molecular formula for a molecule represented as
#' a \code{\link[tidygraph]{tbl_graph}} from \code{tidygraph}.
#'
#' This function calculates the molecular formula of a given
#' \code{\link[tidygraph]{tbl_graph}} object representing a molecular graph. The
#' resulting molecular formula provides information about the types and counts
#' of atoms present in the molecule. This function does not add implicit not
#' remove explicit hydrogen atoms.
#' 
#' @param tbl_graph The input tbl_graph object representing the molecular graph.
#'
#' @return Returns the molecular formula as a character string.
#' 
#' @examples
#' \dontrun{
#'   ## NOT RUN
#'   molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
#'   gr <- sdf2tbl_graph(molec)
#'   mf <- tbl_graph2mf(gr)
#' }
#' 
#' @export
#' 
#' @importFrom tidygraph activate
#' @importFrom dplyr select
#' @importFrom magrittr "%>%"
tbl_graph2mf <- function(tbl_graph) {
  subsMF <- tbl_graph %>%
    activate("nodes") %>%
    as.data.frame() %>%
    select("sym") %>%
    table()

  subsMF <- subsMF[order(names(subsMF))]
  subsMF <- paste(rbind(names(subsMF), subsMF), collapse = "")

  return(subsMF)
}
