
data("lalonde", package="cobalt")
psCols <- c("age", "educ", "married", "nodegree")
treatVal <- "treat"
responseVal <- "re78"  
pairDistVal <- c("age", "married", "educ", "nodegree")
exactVal <- c("educ") 
myBalVal <- c("race")
r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
r2s <- c(0.01)


res1 <- distBalMatch(df=lalonde, treatCol= treatVal,myBalCol = myBalVal, rhoExclude=r1s, rhoBalance=r2s, distList = pairDistVal, exactlist = exactVal, propensityCols = psCols, maxUnMatched = 0.1, caliperOption=NULL, 
                                    toleranceOption=1e-1, maxIter=0, rho.max.f = 10)
## 1. Some successful check on function distBalMatch
test_that("distBalMatch test with correct input", {
  expect_equal(0, sum(unlist(lapply(res1, is.na))))
  expect_equal(res1$version, "Basic")
  expect_equal(length(res1$matchList), 13)
  expect_equal(length(res1$fPair), length(res1$fMarginal))
  expect_equal(length(res1$fPair), length(res1$fExclude))
})

## 2. Check that distBalMatch can handle incorrect input 
test_that("distBalMatch test with wrong input", {
  expect_error(distBalMatch(lalonde, NULL, NULL))
  expect_error(distBalMatch(lalonde, "wrongTreat", "wrongBal"))
})
