## Data generation
set.seed(999)
x1 = rnorm(100, 0, 0.5)
x2 = rnorm(100, 0, 0.1)
x3 = rnorm(100, 0, 1)

x4 = rnorm(100, x1, 0.1)

r1ss <- seq(0.1,50, 5)
r2ss <- seq(0.1,50, 5)

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
d1 <- as.matrix(dist(df["x1"]))
d2 <- as.matrix(dist(df["x2"]))

idx <- 1:length(z)
treatedUnits <- idx[z==1]
controlUnits <- idx[z==0]

d1 <- as.matrix(d1[treatedUnits, controlUnits])
d2 <- as.matrix(d2[treatedUnits, controlUnits])


res6 <- twoDistMatch(df = df, treatCol = "z",  
                     dMat1=d1, dType1= "User", dMat2=d2, dType2="User", myBalCol=c("x5"),  rhoExclude=r1ss, rhoDistance=r2ss, propensityCols = c("x1"), pScores = NULL, ignore = c("y"),  maxUnMatched = 0.1, caliperOption=0.25, 
                     toleranceOption=1e-6, maxIter=3, rho.max.f = 10)
## 0. Check the correctness of the result 
test_that("twoDistMatch test with correct input", {
  expect_equal(names(res6), c("numTreat","rhoList", "matchList", "treatmentCol",
                              "covs", "exactCovs", "idMapping", "b.var",
                              "dataTable", "t", "df", "pair_cost1", "pair_cost2", "version", "fDist1",
                              "fExclude", "fDist2"))
  expect_equal(0, sum(unlist(lapply(res6, is.na))))
  expect_equal(res6$version, "Advanced")
  expect_equal(length(res6$matchList), 236)
  expect_equal(length(res6$fDist1), length(res6$fDist2))
  expect_equal(length(res6$fDist1), length(res6$fExclude))
})



## 1. Tabular summary function works
test_that("Some tabular summary functions only works for the basic version", {
  expect_error(compareMatching(res6))
  summary_table <- generateRhoObj(res6)
  expect_equal(dim(summary_table)[1], 236)
  expect_equal(dim(summary_table)[2], 7)
})

## 2. Some graphical summary function would not work 

test_that("Some graphicall summary functions only works for the basic version", {
  expect_error(generateTVGraph(res6))
  expect_error(generatePairdistanceBalanceGraph(res6))
  expect_error(generatePairdistanceGraph(res6))
})

