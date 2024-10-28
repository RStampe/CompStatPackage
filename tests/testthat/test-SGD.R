library(splines)
set.seed(10)
example_beta <- rnorm(10, 1)
example_lambda <- 0.5
example_x <- seq(0, 10, length.out = 101)
example_y <- rbinom(101, 1, 0.5)

# Test for get_inner_knots function
test_that("get_inner_knots works correctly", {
  result <- CompStatPackage:::get_inner_knots(example_x, length_beta = length(example_beta))

  # Check the length of the inner knots
  expect_equal(length(result), length(example_beta) - 2)


  expect_equal(result, c(0, 1.42857142857143, 2.85714285714286, 4.28571428571429, 5.71428571428571,
                         7.14285714285714, 8.57142857142857, 10))
})

# Test for get_knots function
test_that("get_knots works correctly", {

  example_inner_knots <- seq(0, 10, length.out = 6)
  example_knots <- c(0, 0, 0, example_inner_knots, 10, 10, 10)

  result <- CompStatPackage:::get_knots(inner_knots = example_inner_knots)
  expect_equal(result, example_knots)
})


# Test for calculate_phi function
test_that("calculate_phi works correctly", {

  inner_knots <- CompStatPackage:::get_inner_knots(example_x, length(example_beta))
  knots <- CompStatPackage:::get_knots(inner_knots)
  result <- crossprod(t(CompStatPackage:::calculate_phi(X = example_x, knots = knots)),example_beta)

  # Check result length and type
  expect_equal(length(result), length(example_x))
  expect_true(is.numeric(result))
})



# Test for calculate_Omega function
test_that("calculate_Omega works correctly", {


  inner_knots <- CompStatPackage:::get_inner_knots(example_x, length(example_beta))
  result <- CompStatPackage:::calculate_Omega(inner_knots = inner_knots)

  # Check the result is a matrix or expected structure
  expect_true(is.matrix(result) || is.numeric(result))
  # Additional checks (dimensions, values, etc.)
  expect_equal(ncol(result), length(inner_knots) + 2)
})

# Test for calculate_H function
test_that("calculate_H works correctly", {

  des <- CompStatPackage:::get_design_B_splines(seq(1,2), 5)
  H_object <- CompStatPackage:::function_H_factory(des, y = rep(1,2))

  ## When lambda = 0
  result1 <- H_object$objective(seq(0.1,0.5, length.out = 5), lambda = 0)
  expect_equal(as.numeric(result1), 0.5592368, tolerance = 1e-6)

  result2 <- H_object$objective(seq(0.1,0.5, length.out = 5), lambda = 1)
  expect_true(result1 < result2)
})





# Test for calculate_H function
test_that("deriv of H works correctly", {

  lambda = 1e-2
  X = rnorm(20)

  y = rnorm(20)
  p = 10
  length_beta= 10
  beta = rnorm(p)


  des <- CompStatPackage:::get_design_B_splines(X, p)
  H_object <- CompStatPackage:::function_H_factory(des, y)

  error_margin <- range(H_object$gradient(beta, lambda) - grad(function(beta) H_object$objective(beta, lambda), beta))

  expect_equal(error_margin[1], 0, tolerance = 1e-9)
  expect_equal(error_margin[2], 0, tolerance = 1e-9)
})




# Test for calculate_H function
test_that("Hessian_numeric of H works correctly", {

  lambda = 1e-2
  X = rnorm(20)

  y = rnorm(20)
  p = 10
  length_beta= 10
  beta = rnorm(p)


  des <- CompStatPackage:::get_design_B_splines(X, p)
  H_object <- CompStatPackage:::function_H_factory(des, y)

  error_margin <- range(H_object$Hessian(beta, lambda) - hessian(function(beta) H_object$objective(beta, lambda), beta))

  expect_equal(error_margin[1], 0, tolerance = 1e-8)
  expect_equal(error_margin[2], 0, tolerance = 1e-8)
})

