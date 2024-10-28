library(testthat)
library(numDeriv)

# Define parameters for the tests
parameters <- c(mu = 0, sigma2 = 1, nu = 3)
x <- c(-1, 0, 1)

# Create EM setup
em_setup <- get_em_setup(x)
test_that("e_step function returns expected values", {
  e_step_values <- em_setup$e_step(parameters)
  expect_type(e_step_values, "double")
  expect_equal(length(e_step_values), length(x))
  expect_true(all(e_step_values > 0))
})

test_that("m_step function computes updated parameters correctly", {
  E <- c(1.2, 0.8, 1.1)
  updated_params <- em_setup$m_step(E, parameters)
  expect_type(updated_params, "double")
  expect_named(updated_params, c("mu", "sigma2", "nu"))
  expect_true(updated_params["sigma2"] > 0)  # Variance should be positive
})


test_that("loglike_grad matches numeric gradient", {
  set.seed(100)
  data <- sample_data_for_em(1000)
  fct <- get_gd_setup(data$x)
  par_test <- c(mu = 1, sigma2 = 4, nu = 2)
  grad_analytic <- fct$H_grad(par_test)
  grad_numeric <- grad(function(par) fct$H(par), par_test)

  expect_equal(unname(grad_analytic)[1:2], unname(grad_numeric)[1:2], tolerance = 1e-6)
})





