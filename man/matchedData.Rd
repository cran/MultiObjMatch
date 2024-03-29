% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiObjMatch.R
\name{matchedData}
\alias{matchedData}
\title{Get matched dataframe}
\usage{
matchedData(matchingResult, match_num)
}
\arguments{
\item{matchingResult}{an object returned by the main matching function
distBalMatch}

\item{match_num}{Integer index of match that the user want to extract paired
observations from}
}
\value{
dataframe that contains only matched pair data
}
\description{
A function that returns the dataframe that contains only matched
pairs from the original data frame with specified match index
}
\examples{
## Generate Matches
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
matchedData(matchResult1, 1)
}
