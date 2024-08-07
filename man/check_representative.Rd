% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiObjMatch.R
\name{check_representative}
\alias{check_representative}
\title{Check the representativeness of matched treated units}
\usage{
check_representative(matching_result, match_num = NULL)
}
\arguments{
\item{matching_result}{the matching result returned by either dist_bal_match
or two_dist_match.}

\item{match_num}{(optional) Integer index of match that the user want to extract paired
observations from. NULL by default, which will generate a table for all the matches.}
}
\value{
a summary table of SMDs of the key covariates between the whole
treated units and the matched treated units.
}
\description{
Summary function to compare SMD of the key covariates in matched
and the full set of treated units.
}
