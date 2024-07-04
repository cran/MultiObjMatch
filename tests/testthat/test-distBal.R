
data("lalonde", package="cobalt")
psCols <- c("age", "educ", "married", "nodegree")
treatVal <- "treat"
responseVal <- "re78"  
pairDistVal <- c("age", "married", "educ", "nodegree")
exactVal <- c("educ") 
myBalVal <- c("race")
r1s <- c( 0.1, 0.3, 0.5, 0.7, 0.9,1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7)
r2s <- c(0.01)


res1 <- dist_bal_match(data=lalonde, treat_col= treatVal, 
                     marg_bal_col = myBalVal, 
                     exclusion_penalty=r1s, 
                     balance_penalty=r2s, 
                     dist_col = pairDistVal, 
                     exact_col = exactVal, 
                     propensity_col = psCols, 
                     max_unmatched = 0.1, 
                     caliper_option=NULL, 
                     tol=1e-1, max_iter=0, rho_max_factor = 10)
## 1. Some successful check on function dist_bal_match
test_that("dist_bal_match test with correct input", {
  expect_equal(0, sum(unlist(lapply(res1, is.na))))
  expect_equal(res1$version, "Basic")
  expect_equal(length(res1$matchList), 13)
  expect_equal(length(res1$fPair), length(res1$fMarginal))
  expect_equal(length(res1$fPair), length(res1$fExclude))
})

## 2. Check that dist_bal_match can handle incorrect input 
test_that("dist_bal_match test with wrong input", {
  expect_error(dist_bal_match(lalonde, NULL, NULL))
  expect_error(dist_bal_match(lalonde, "wrongTreat", "wrongBal"))
})
