% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiObjMatch.R
\name{compareMatching}
\alias{compareMatching}
\title{Generate covariate balance in different matches}
\usage{
compareMatching(
  matchingResult,
  covList = NULL,
  display.all = TRUE,
  stat = "mean.diff"
)
}
\arguments{
\item{matchingResult}{an object returned by the main matching function
distBalMatch}

\item{covList}{(optional) factor of names of covariates that we want to
evaluate covariate balance on; default is NULL. When set to NULL, the
program will compare the covariates that have been used to construct a
propensity model.}

\item{display.all}{(optional) boolean value of whether to display all the
matches; default is TRUE, where matches at each quantile is displayed}

\item{stat}{(optional) character of the name of the statistic used for
measuring covariate balance; default is "mean.diff". This argument is the
same as used in "cobalt" package, see: \link[cobalt]{bal.tab}}
}
\value{
a dataframe that shows covariate balance in different matches
}
\description{
This is a wrapper function for use in evaluating covariate
balance across different matches. The function calls \code{compareTables} on the
output from the function \code{generateBalanceTable}. It only works for 'Basic'
version of matching (using \code{distBalMatch}).
}
\examples{
## Generate matches 
data("lalonde", package="cobalt")
psCols <- c("age", "educ", "married", "nodegree")
treatVal <- "treat"
responseVal <- "re78"
pairDistVal <- c("age", "married", "educ", "nodegree")
exactVal <- c("educ")
myBalVal <- c("race")
r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
r2s <- c(0.01)
matchResult <- distBalMatch(df=lalonde, treatCol=treatVal, myBalCol=myBalVal,
rhoExclude =r1s, rhoBalance=r2s,
distList=pairDistVal, exactlist=exactVal,
propensityCols = psCols,ignore = c(responseVal), maxUnMatched = 0.1,
caliperOption=NULL, toleranceOption=1e-1, maxIter=0, rho.max.f = 10)

## Generate table for comparing matches
compareMatching(matchResult, display.all = TRUE)
}
