## This is the file for testing the helper functions in the package MultiObjMatch

data("lalonde", package="cobalt")
psCols <- c("age", "educ", "married", "nodegree")
treatVal <- "treat"
responseVal <- "re78"  
pairDistVal <- c("age", "married", "educ", "nodegree")
exactVal <- c("educ") 
myBalVal <- c("race")
r1s <- seq(0.1, 10, 0.5)
r2s <- c(0.001)


res1 <- distBalMatch(df=lalonde, treatCol= treatVal,myBalCol = myBalVal, rhoExclude=r1s, rhoBalance=r2s, distList = pairDistVal, exactlist = exactVal, propensityCols = psCols, ignore = c(responseVal), maxUnMatched = 0.1, caliperOption=NULL, 
                     toleranceOption=1e-1, maxIter=0, rho.max.f = 10)

## 1. Tabular summary function works
test_that("Tabular summary function test", {
  matching_comparison <- compareMatching(res1)
  matching_comparison_full <- compareMatching(res1, display.all = FALSE)
  matching_comparison_partial <- compareMatching(res1, covList = c("age", "married"))
  expect_equal(class(matching_comparison), "data.frame")
  expect_equal(class(matching_comparison_full), "data.frame")
  expect_equal(class(matching_comparison_partial), "data.frame")
  expect_equal(dim(matching_comparison)[2], dim(matching_comparison_partial)[2])
  expect_equal(dim(matching_comparison)[1], dim(matching_comparison_full)[1])
  expect_equal(length(c("age", "married")), dim(matching_comparison_partial)[1])
})

## 2. Test the function that generate the objective function values and penalty
test_that("Objective function summary", {
  summary_table <- generateRhoObj(res1)
  expect_equal(dim(summary_table)[1], 20)
  expect_equal(dim(summary_table)[2], 7)
})

## 3. Test graphic helper function 
test_that("Graphic helper function", {
  
  expect_error(visualize(res1, "f1"))
  expect_error(visualize(res1, x_axis = "f1", y_axis = "f2"))
  
})





