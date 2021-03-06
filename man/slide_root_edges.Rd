% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{slide_root_edges}
\alias{slide_root_edges}
\title{Repartitions lengths along edges that descend from root node so
that they are equal. Useful after rooting operations that result
in branches with 0 length}
\usage{
slide_root_edges(phy)
}
\arguments{
\item{phy}{A phylo object}
}
\value{
A phylo object with modified edge lengths
}
\description{
Repartitions lengths along edges that descend from root node so
that they are equal. Useful after rooting operations that result
in branches with 0 length
}
