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


res1 <- dist_bal_match(data=lalonde, treat_col= treatVal,marg_bal_col = myBalVal, 
                     exclusion_penalty=r1s, balance_penalty=r2s, 
                     dist_col = pairDistVal, 
                     exact_col= exactVal, 
                     propensity_col = psCols, 
                     ignore_col = c(responseVal), 
                     max_unmatched = 0.1, caliper_option=NULL, 
                     tol=1e-1, max_iter=0, rho_max_factor = 10)

## 1. Tabular summary function works
test_that("Tabular summary function test", {
  matching_comparison <- compare_matching(res1)
  matching_comparison_full <- compare_matching(res1, display_all = FALSE)
  matching_comparison_partial <- compare_matching(res1, cov_list = c("age", "married"))
  expect_equal(class(matching_comparison), "data.frame")
  expect_equal(class(matching_comparison_full), "data.frame")
  expect_equal(class(matching_comparison_partial), "data.frame")
  expect_equal(dim(matching_comparison)[2], dim(matching_comparison_partial)[2])
  expect_equal(dim(matching_comparison)[1], dim(matching_comparison_full)[1])
  expect_equal(length(c("age", "married")), dim(matching_comparison_partial)[1])
})

## 2. Test the function that generate the objective function values and penalty
test_that("Objective function summary", {
  summary_table <- get_rho_obj(res1)
  expect_equal(dim(summary_table)[1], 20)
  expect_equal(dim(summary_table)[2], 7)
})

## 3. Test graphic helper function 
test_that("Graphic helper function", {
  
  expect_error(visualize(res1, "f1"))
  expect_error(visualize(res1, x_axis = "f1", y_axis = "f2"))
  
})





