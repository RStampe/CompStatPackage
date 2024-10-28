library(testthat)
library(Rcpp)

test_that("epanechnikov_kernel returns correct values", {
  expect_equal(CompStatPackage:::epanechnikov_kernel(c(0,1,10,0.5)), c(0.75,0,0,0.5625))
  expect_equal(CompStatPackage:::epanechnikov_kernel(c(0,1,10,0.5)), c(0.75,0,0,0.5625))
})

test_that("gaussian_kernel returns correct values", {
  expect_equal(CompStatPackage:::gaussian_kernel(0), 1 / sqrt(2 * pi))
  expect_equal(CompStatPackage:::gaussian_kernel(1), (1 / sqrt(2 * pi)) * exp(-0.5))
})


# test_that("kernel evaluators returns correct density estimates", {
#   grid <- c(0, 1, 2, 3)
#   vec_obs <- c(0.5, 1.5, 2, 1)
#   result <- CompStatPackage:::kernel_vec(grid, vec_obs, bandwidth = 0.5, CompStatPackage:::epanechnikov_kernel_cpp)
#   result_cpp <- CompStatPackage:::kernel_cpp_epanechnikov(grid, vec_obs, bandwidth = 0.5)
#
#   expect_type(result, "list")
#   expect_equal(result$x, grid)
#   expect_equal(result$y, c(0.000, 0.375, 0.375, 0.000))
#   expect_equal(result, result_cpp)
# })


test_that("kernel_density_estimation computes density estimates", {
  vec_obs <- c(0.5, 1.5, 2.5)
  result <- CompStatPackage:::kernel_density_estimation(vec_obs, bandwidth = 1, kernel_calculator = CompStatPackage:::kernel_vec, kernel_function = CompStatPackage:::gaussian_kernel)

  expect_s3_class(result, "density")
  expect_equal(result$n, length(vec_obs))
  expect_true(all(result$y >= 0))
})


test_that("bandwidth_silvermans_rule calculates bandwidth correctly", {
  vec_obs <- c(1, 2, 3, 4, 5)
  bw <- CompStatPackage:::bandwidth_silvermans_rule(vec_obs)

  expect_gt(bw, 0)
})

test_that("get bw cv is correct", {
  set.seed(10)
  n <- 1000
  norm_data <- rnorm(n)

  bw <- CompStatPackage:::get_cv_bw(norm_data)

  expect_equal(bw, 0.58, tolerance = 0.02)
})

test_that("get plug-in is correct", {
  set.seed(10)
  n <- 1000
  norm_data <- rnorm(n)

  bw <- CompStatPackage:::calculate_oracle_bw_epanechnikov(vec_obs = norm_data)

  expect_equal(bw, 0.3, tolerance = 0.5)

})

test_that("get bw silverman is correct", {
  set.seed(10)
  n <- 1000
  norm_data <- rnorm(n)

  bw <- CompStatPackage:::bandwidth_silvermans_rule(norm_data)

  expect_equal(bw, 0.222396, tolerance = 0.0001)

})

test_that("get bw ucv is correct", {
  set.seed(10)
  n <- 1000
  norm_data <- rnorm(n)

  bw <- CompStatPackage:::get_ucv_bw(norm_data)

  expect_equal(bw, 0.483, tolerance = 0.001)
})

# test_that("binning_gives_correct_results", {
#   set.seed(10)
#   vec = rnorm(10000)
#   n <- 10
#   expected <- c(7e-04, 0.0107, 0.0583, 0.1794, 0.3071, 0.2704, 0.1369, 0.0308,
#                 0.0052, 5e-04)
#
#   bin_Vals_cpp <- CompStatPackage:::binning_cpp(vec, lower = min(vec), upper = max(vec), number_of_bins = n)
#   expect_equal(bin_Vals_cpp, expected)
# })

