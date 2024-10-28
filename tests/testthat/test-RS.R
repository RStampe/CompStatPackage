library(testthat)

test_that("targets calculates correctly", {
  grid <- seq(-100,200, length.out = 10000)
  import_poisson_data()
  vectorize = log_target_density_vectorize(grid)
  cpp = log_target_density_cpp(grid)


  expect_equal(vectorize, cpp)
})

test_that("get_proposal functions return correct value", {
  import_poisson_data()
  setup <- get_gaussian_envelope(log_target_density_vectorize, type = "t")
  n <- 100

  vectorize <- get_proposals.vectorize(setup, n)


  expect_equal(length(vectorize), 100)
  expect_equal(mean(vectorize), 0.24, tolerance = 0.2)


  cpp <- get_proposals.cpp(setup, n)


  expect_equal(length(cpp), 100)
  expect_equal(mean(cpp), 0.24, tolerance = 0.2)
})

test_that("get_accepts functions return correct value",{
  import_poisson_data()
  set.seed(18123)
  setup <- get_gaussian_envelope(log_target_density_vectorize, type = "t")
  proposals <- get_proposals.vectorize(setup, 100)

  vectorize <- get_accepts.vectorize(setup, proposals)

  expect_equal(length(vectorize), 100)
  expect_equal(mean(vectorize), 0.93, tolerance = 0.05)


  cpp <- get_accepts.cpp(setup, proposals)

  expect_equal(length(cpp), 100)
  expect_equal(mean(cpp), 0.93, tolerance = 0.05)
})
