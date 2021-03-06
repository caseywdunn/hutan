% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{distance_from_tip}
\alias{distance_from_tip}
\title{For each node in the tree (including internal nodes and tips)
get the shortest distance to a descendant node. Values for tip
nodes should be 0.}
\usage{
distance_from_tip(phy)
}
\arguments{
\item{phy}{A phylo object}
}
\value{
A vector with elements corresponding to each node
}
\description{
For each node in the tree (including internal nodes and tips)
get the shortest distance to a descendant node. Values for tip
nodes should be 0.
}
