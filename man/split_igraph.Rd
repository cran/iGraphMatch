% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_igraph.R
\name{split_igraph}
\alias{split_igraph}
\title{Split an igraph object into aligned graphs by attribute}
\usage{
split_igraph(g, e_attr, strip_vertex_attr = FALSE)
}
\arguments{
\item{g}{An igraph object}

\item{e_attr}{the name of an edge attribute in g}

\item{strip_vertex_attr}{Whether to remove all vertex
attribute from the new graphs}
}
\value{
A named list of igraph objects with names corresponding to the values of
 the edge attributes.
}
\description{
Given an igraph object and an edge attribute, this function
finds all unique values of the edge attribute in the graph
and returns a list of igraph objects on the same vertex set
where each element of the list has a graph containing only
those edges with specified attributed.
}
\examples{
g <- igraph::sample_gnm(20, 60)
igraph::E(g)$color <-
  sample(c("red", "green"), 60, replace = TRUE)
split_igraph(g, "color")

}
