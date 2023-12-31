% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GraphProcessing.R
\name{split_graph_at}
\alias{split_graph_at}
\title{Split molecular graph by removing a specified atom type}
\usage{
split_graph_at(tbl_graph, splitpos = "P")
}
\arguments{
\item{tbl_graph}{The input tbl_graph object representing the molecular graph.}

\item{splitpos}{A character specifying the atom type to be removed by the
splitting process.}
}
\value{
Returns a list of tbl_graph objects, each representing a separate
  molecular graph after the split.
}
\description{
This function splits the given \code{\link[tidygraph]{tbl_graph}} object by
removing all by-type specified atoms from the graph. The resulting graph(s)
consist(s) of (multiple) connected component(s), each representing a separate
molecular graph after the split.
}
\examples{
\dontrun{
  ## NOT RUN
  molec <- ChemmineR::smiles2sdf("CCP(CO)CF(F)(F)")
  gr <- sdf2tbl_graph(molec)
  subs <- split_graph_at(gr, splitpos = "P")
}

}
