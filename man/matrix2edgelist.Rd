% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tradeoff_functions.R
\name{matrix2edgelist}
\alias{matrix2edgelist}
\title{Helper function to convert matrix to list}
\usage{
matrix2edgelist(mat)
}
\arguments{
\item{mat}{matrix representation of distances between treated and control
units}
}
\value{
list of vector representation of distances
}
\description{
Convert between a matrix representation of distances between
treated and control units and a list of vectors (default format for
build.dist.struct function in rcbalance package)
}
