% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tradeoff_functions.R
\name{callrelax}
\alias{callrelax}
\title{Call relax on the network}
\usage{
callrelax(net, solver = "rlemon")
}
\arguments{
\item{net}{the network structure}

\item{solver}{(optional) the solver; by default, "rlemon"}
}
\value{
list of the result from the call to relax solver
}
\description{
this function is copied from the rcbalance package
}
