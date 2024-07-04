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


res6 <- two_dist_match(data = df, treat_col = "z",  
                     dist1_matrix=d1, dist1_type= "User", 
                     dist2_matrix=d2, dist2_type="User", 
                     marg_bal_col=c("x5"),  
                     exclusion_penalty=r1ss, 
                     dist2_penalty=r2ss, 
                     propensity_col = c("x1"), 
                     pscore_name = NULL, 
                     ignore_col = c("y"),  
                     max_unmatched = 0.1, 
                     caliper_option=0.25, 
                     tol=1e-6, 
                     max_iter=0, 
                     rho_max_factor = 10)
## 0. Check the correctness of the result 
test_that("two_dist_match test with correct input", {
  expect_equal(names(res6), c("numTreat","rhoList", "matchList", "treatmentCol",
                              "covs", "exactCovs", "idMapping", "b.var",
                              "dataTable", "t", "df", "pair_cost1", "pair_cost2", "version", "fDist1",
                              "fExclude", "fDist2"))
  expect_equal(0, sum(unlist(lapply(res6, is.na))))
  expect_equal(res6$version, "Advanced")
  expect_equal(length(res6$matchList), 100)
  expect_equal(length(res6$fDist1), length(res6$fDist2))
  expect_equal(length(res6$fDist1), length(res6$fExclude))
})



## 1. Tabular summary function works
test_that("Some tabular summary functions only works for the basic version", {
  expect_error(compare_matching(res6))
  summary_table <- get_rho_obj(res6)
  expect_equal(dim(summary_table)[1], 100)
  expect_equal(dim(summary_table)[2], 7)
})

## 2. Some graphical summary function would not work 

test_that("Some graphicall summary functions only works for the basic version", {
  expect_error(get_tv_graph(res6))
  expect_error(get_pairdist_graph(res6))
  expect_error(get_pairdist_balance_graph(res6))
})

