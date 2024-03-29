% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/largest_common_cc.R
\name{largest_common_cc}
\alias{largest_common_cc}
\alias{largest_cc}
\title{Find the largest common connected subgraph (LCCS) of two graphs}
\usage{
largest_common_cc(A, B, min_degree = 1)

largest_cc(A)
}
\arguments{
\item{A}{A matrix or an igraph object. See \link{check_graph}. Must be single-layer.}

\item{B}{A matrix or an igraph object. See \link{check_graph}. Must be single-layer.}

\item{min_degree}{A number. Defines the level of connectedness of the
obtained largest common connected subgraph. The induced subgraph is
a graph with a minimum vertex-degree of at least min_degree.}
}
\value{
\code{largest_common_cc} returns the common largest connected subgraphs of
two aligned graphs in the igraph object form and a logical vector indicating which vertices in
the original graphs remain in the induced subgraph.
}
\description{
Find the largest common connected subgraphs of
two matched graphs, which is an induced connected subgraph of both graphs
that has as many vertices as possible.
The \code{largest_cc} function returns the largest connected subgraph of a single graph.
}
\examples{
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.7, p =  0.2)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
# put no constraint on the minimum degree of the common largest conncect subgraph
lccs1 <- largest_common_cc(g1, g2, min_degree = 1)
# induced subgraph
lccs1$g1
lccs1$g2
# label of vertices of the induced subgraph in the original graph
igraph::V(g1)[lccs1$keep]

# obtain a common largest connect subgraph with each vertex having a minimum degree of 3
lccs3 <- largest_common_cc(g1, g2, min_degree = 3)

g <- igraph::sample_gnp(100, .01)
lcc <- largest_cc(g)
# induced subgraph
lcc$g
# label of vertices of the induced subgraph in the original graph
igraph::V(g)[lcc$keep]

}
