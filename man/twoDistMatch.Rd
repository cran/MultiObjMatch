% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiObjMatch.R
\name{twoDistMatch}
\alias{twoDistMatch}
\title{Optimal tradeoffs among two distances and exclusion}
\usage{
twoDistMatch(
  dType1 = "user",
  dType2 = "user",
  dMat1 = NULL,
  df = NULL,
  dMat2 = NULL,
  treatCol = NULL,
  distList1 = NULL,
  distList2 = NULL,
  rhoExclude = c(1),
  rhoDistance = c(1, 2, 3),
  myBalCol = NULL,
  exactlist = NULL,
  propensityCols = NULL,
  pScores = NULL,
  ignore = NULL,
  maxUnMatched = 0.25,
  caliperOption = NULL,
  toleranceOption = 0.01,
  maxIter = 0,
  rho.max.f = 10
)
}
\arguments{
\item{dType1}{One of ("Euclidean", "Mahalanobis", "user") indicating the type
of distance that are used for the first distance objective functions. NULL
by default.}

\item{dType2}{One of ("Euclidean", "Mahalanobis", "user")  charactor
indicating the type of distance that are used for the second distance
objective functions. NULL by default.}

\item{dMat1}{(optional) matrix object that represents the distance matrix
using the first distance measure; \code{dType} must be passed in as "user" if
dMat is non-empty}

\item{df}{(optional) data frame that contain columns indicating treatment,
outcome and covariates}

\item{dMat2}{(optional) matrix object that represents the distance matrix
using the second distance measure; \code{dType1} must be passed in as "user" if
dMat is non-empty}

\item{treatCol}{(optional) character, name of the column indicating treatment
assignment.}

\item{distList1}{(optional) character vector names of the variables used for
calculating covariate distance using first distance measure specified by
dType}

\item{distList2}{(optional) character vector, names of the variables used for
calculating covariate distance using second distance measure specified by
dType1}

\item{rhoExclude}{(optional) numeric vector, penalty values associated with
the distance specified by \code{dMat} or \code{dType}. Default value is c(1).}

\item{rhoDistance}{(optional) numeric vector, penalty values associated with
the distance specified by \code{dMat1} or \code{dType1}. Default value is c(1,2,3).}

\item{myBalCol}{(optional) character, column name of the variable on which to
evaluate balance.}

\item{exactlist}{(optional) character vector, names of the variables on which
to match exactly; NULL by default.}

\item{propensityCols}{character vector, names of columns on which to fit a
propensity score model.}

\item{pScores}{(optional) character, name of the column containing fitted
propensity scores; default is NULL.}

\item{ignore}{(optional) character vector of variable names that should be
ignored when constructing the internal matching. NULL by default.}

\item{maxUnMatched}{(optional) numeric, maximum proportion of unmatched units
that can be accepted; default is 0.25.}

\item{caliperOption}{(optional) numeric, the propensity score caliper value
in standard deviations of the estimated propensity scores; default is NULL,
which is no caliper.}

\item{toleranceOption}{(optional) numeric, tolerance of close match distance;
default is 1e-2.}

\item{maxIter}{(optional) integer,  maximum number of iterations to use in
searching for penalty combintions that improve the matching; default is 0.}

\item{rho.max.f}{(optional) numeric, the scaling factor used in proposal for
rhos; default is 10.}
}
\value{
a named list whose elements are:
\itemize{
\item "rhoList": list of penalty
combinations for each match
\item "matchList": list of matches indexed by
number
\item "treatmentCol": character of treatment variable
\item "covs":character vector of names of the variables used for calculating within-pair
distance
\item "exactCovs": character vector of names of variables that we want
exact or close match on
\item "idMapping": numeric vector of row indices for
each observation in the sorted data frame for internal use
\item "stats": data
frame of important statistics (total variation distance) for variable on
which marginal balance is measured
\item "b.var": character, name of variable
on which marginal balance is measured (left NULL since no balance
constraint is imposed here).
\item "dataTable": data frame sorted by treatment
value
\item "t": a treatment vector
\item "df": the original dataframe input by the
user
\item "pair_cost1": list of pair-wise distance sum using the first
distance measure
\item "pair_cost2": list of pair-wise distance sum using the
second distance measure
\item "version": (for internal use) the version of the
matching function called; "Basic" indicates the matching comes from
distBalMatch and "Advanced" from twoDistMatch.
\item "fDist1": a vector of
values for the first objective function; it corresponds to the pair-wise
distance sum according to the first distance measure.
\item "fExclude": a
vector of values for the second objective function; it corresponds to the
number of treated units being unmatched.
\item "fDist2": a vector of values for
the third objective function; it corresponds to the pair-wise distance sum
corresponds to the
}
}
\description{
Explores tradeoffs among three objective functions in
multivariate matching: sums of two different user-specified  covariate
distances within matched pairs, and the number of treated units included in
the match.
}
\details{
Matched designs generated by this function are Pareto optimal for
the three objective functions.  The degree of relative emphasis among the
three objectives in any specific solution is controlled by the penalties,
denoted by Greek letter rho. Larger values for the penalties associated with
the two distances correspond to increased emphasis close matching on these
distances, at the possible cost of excluding more treated units. Additional
details:
\itemize{
\item Users may either specify their own distance matrices (specifying
the \code{user} option in \code{dType1} and/or \code{dType2} and supplying arguments to
\code{dMat1} and/or \code{dMat2} respectively) or ask the function to create
Mahalanobis or Euclidean distances on sets of covariates specified by the
\code{distList1} and \code{distList2} arguments. If \code{dType1} or \code{dType2} is not
specified, if one of these is set to \code{user} and the corresponding \code{dMat1}
argument is not provided, or if one is NOT set to \code{user} and the
corresponding \code{distList1} argument is not provided, an error will result.
\item User-specified distance matrices passed to \code{dMat1} or \code{dMat2} should have row
count equal to the number of treated units and column count equal to the
number of controls.
\item If the \code{caliperOption} argument is specified, a
propensity score caliper will be imposed, forbidding matches between units
more than a fixed distance apart on the propensity score.  The caliper will
be based either on a user-fit propensity score, identified in the input
dataframe by argument \code{pScores}, or by an internally-fit propensity score
based on logistic regression against the variables named in \code{psoreCols}.  If
\code{caliperOption} is non-NULL and neither of the other arguments is specified
an error will result.
\item \code{toleranceOption} controls the precision at which the
objective functions is evaluated. When matching problems are especially large
or complex it may be necessary to increase toleranceOption in order to
prevent integer overflows in the underlying network flow solver; generally
this will be suggested in appropariate warning messages.
\item While by default
tradeoffs are only assessed at penalty combinations provided by the user, the
user may ask for the algorithm to search over additional penalty values in
order to identify additional Pareto optimal solutions. \code{rho.max.f} is a
multiplier applied to initial penalty values to discover new solutions, and
setting it larger leads to wider exploration; similarly, \code{maxIter} controls
how long the exploration routine runs, with larger values leading to more
exploration.
}
}
\examples{
x1 = rnorm(100, 0, 0.5)
x2 = rnorm(100, 0, 0.1)
x3 = rnorm(100, 0, 1)
x4 = rnorm(100, x1, 0.1)
r1ss <- seq(0.1,50, 10)
r2ss <- seq(0.1,50, 10)
x = cbind(x1, x2, x3,x4)
z = sample(c(rep(1, 50), rep(0, 50)))
e1 = rnorm(100, 0, 1.5)
e0 = rnorm(100, 0, 1.5)
y1impute = x1^2 + 0.6*x2^2 + 1 + e1
y0impute = x1^2 + 0.6*x2^2 + e0
treat = (z==1)
y = ifelse(treat, y1impute, y0impute)
names(x) <- c("x1", "x2", "x3", "x4")
df <- data.frame(cbind(z, y, x))
df$x5 <- 1
names(x) <- c("x1", "x2", "x3", "x4")
df <- data.frame(cbind(z, y, x))
df$x5 <- 1
d1 <- as.matrix(dist(df["x1"]))
d2 <- as.matrix(dist(df["x2"]))
idx <- 1:length(z)
treatedUnits <- idx[z==1]
controlUnits <- idx[z==0]
d1 <- as.matrix(d1[treatedUnits, controlUnits])
d2 <- as.matrix(d2[treatedUnits, controlUnits])
matchResult1 <- twoDistMatch(df, "z", "y", dMat1=d1, dType1= "User", dMat2=d2,
dType2="User", myBalCol=c("x5"), rhoExclude=r1ss, rhoDistance=r2ss,
propensityCols = c("x1")) 
}
\seealso{
Other main matching function: 
\code{\link{distBalMatch}()}
}
\concept{main matching function}
