#' Generate a sample of data with parameters
#'
#' This function generates a random sample of data `x` from a normal distribution
#' with mean `mu`, variance `sigma2`, and degrees of freedom `nu`. It also provides
#' a `w` vector of chi-squared distributed values.
#'
#' @param n An integer specifying the number of samples to generate.
#' @param parameters A numeric vector of length three for `mu`, `sigma2`, and `nu`.
#' @return A list with elements:
#'   \item{x}{A numeric vector of sampled data.}
#'   \item{w}{A numeric vector of chi-squared distributed values.}
#' @export
sample_data_for_em <- function(n, parameters = c(0.5, 1, 1)) {
  mu <- parameters[1]
  sigma2 <- parameters[2]
  nu <- parameters[3]

  w <- rchisq(n, df = nu)
  x <- rnorm(n, mu, sqrt(nu * sigma2 / w))
  list(x = x, w = w)
}

#' Optimization using EM or GD
#'
#' Provides a unified interface for optimizing a setup using either the expectation-maximization (EM)
#' or gradient descent (GD) methods.
#'
#' @param x A numeric vector of data to optimize over.
#' @param type A character string specifying the optimization type, either "em" for EM or "gd" for GD.
#' @param configuration A list of optimization configuration settings, including `initial_parameters`,
#'   `stopping_tolerance`, and `gamma`.
#' @param do_trace A logical, if `TRUE`, enables tracing during optimization.
#' @return An object of class `optim_result` containing optimized parameters and configuration details.
#' @export
perform_em_or_gd <- function(x, type = "em", configuration = list(), do_trace = FALSE) {
  default_configuration <- list(
    initial_parameters = c(mu = 1, sigma2 = 1, nu = 5),
    stopping_tolerance = 1e-16,
    gamma = 0.01
  )

  configuration <- modifyList(default_configuration, configuration)


  if (exists("my_debug")) if (my_debug) browser()
  type <- match.arg(type, c("em", "gd"))
  setup <- if (type == "em") {
    get_em_setup(x)
  } else {
    get_gd_setup(x)
  }

  structure(optimize_this(setup, configuration, do_trace), class = "optim_result")
}

#' Generic Optimize Function
#'
#' Dispatches `optimize_this` to the appropriate method based on the object class.
#'
#' @param setup An object representing the setup structure for optimization.
#' @param con A list of optimization configurations.
#' @param do_trace A logical, if `TRUE`, enables tracing during optimization.
#' @param ... Additional arguments passed to methods.
optimize_this <- function(setup, con, do_trace, ...) {
  UseMethod("optimize_this")
}

#' Set Up EM Optimization
#'
#' Constructs the expectation-maximization (EM) optimization setup.
#'
#' @param x A numeric vector of observed data.
#' @return A list structure for EM containing `e_step`, `m_step`, and `x`.
get_em_setup <- function(x) {
  e_step <- function(parameters) {
    mu <- parameters[1]
    sigma2 <- parameters[2]
    nu <- parameters[3]

    numerator <- (nu + 1) * nu * sigma2
    denominator <- nu * sigma2 + (x - mu)^2
    numerator / denominator
  }

  m_step <- function(E, parameters) {
    nu <- unname(parameters)[3]

    mu_hat <- sum(x * E) / sum(E)
    sigma2_hat <- mean(E * (x - mu_hat)^2) / nu

    c(mu = mu_hat, sigma2 = sigma2_hat, nu = nu)
  }

  structure(
    list(e_step = e_step, m_step = m_step, x = x),
    class = "em"
  )
}


#' EM Optimization Method
#'
#' Executes the optimization using the expectation-maximization (EM) approach.
#'
#' @param setup An object of class `em` created by `get_em_setup`.
#' @param con A list of optimization configurations.
#' @param do_trace A logical, if `TRUE`, enables tracing during optimization.
#' @return A list with optimized parameters.
optimize_this.em <- function(setup, con, do_trace) {
  epsilon <- con$stopping_tolerance
  parameters <- con$initial_parameters

  epsilon_squared <- epsilon^2
  not_converged <- TRUE
  while (not_converged) {
    old_parameters <- parameters
    E <- setup$e_step(old_parameters)
    parameters <- setup$m_step(E, old_parameters)
    not_converged <- crossprod(parameters - old_parameters) >
      epsilon_squared * (crossprod(old_parameters) + epsilon)^2
  }
  list(
    parameters = parameters,
    con = con,
    setup
  )
}







#' Set Up GD Optimization
#'
#' Constructs the gradient descent (GD) optimization setup.
#'
#' @param x A numeric vector of observed data.
#' @return A list structure for GD containing `H`, `H_grad`, and `x`.
get_gd_setup <- function(x) {
  negative_log_likelihood <- function(x, parameters) {
    mu <- parameters[1]
    sigma2 <- parameters[2]
    nu <- parameters[3]
    n <- length(x)


    -length(x) / 2 * log(sigma2) -
      (nu + 1) / 2 * sum(log(1 + (x - mu)^2 / (nu * sigma2)))
    #  -n/2*log(sigma2) - (nu + 1)/2* sum(log(1 + (x - mu)^2 /(nu * sigma2)))
  }

  negative_log_likelihood_gradient <- function(x, parameters) {
    mu <- parameters[1]
    sigma2 <- parameters[2]
    nu <- parameters[3]
    n <- length(x)

    # d_mu <- -(nu+1)* sum((x-mu) /((x-mu)^2+nu* sigma2))
    # d_sigma2 <- -nu/2*sum(((x-mu)^2-sigma2)/(sigma2*(nu*sigma2+ (x-mu)^2)))
    #
    d_mu <- (nu + 1) * sum((x - mu) / ((x - mu)^2 + nu * sigma2))
    d_sigma2 <- nu / 2 * sum(((x - mu)^2 - sigma2) /
                               (sigma2 * (nu * sigma2 + (x - mu)^2)))


    c(d_mu = d_mu, d_sigma2 = d_sigma2, d_nu = 0)
  }

  structure(
    list(
      H = function(parameters) negative_log_likelihood(x, parameters),
      H_grad = function(parameters) negative_log_likelihood_gradient(x, parameters),
      x = x
    ),
    class = "gd"
  )
}

#' GD Optimization Method
#'
#' Executes the optimization using the gradient descent (GD) approach.
#'
#' @param setup An object of class `gd` created by `get_gd_setup`.
#' @param con A list of optimization configurations.
#' @param do_trace A logical, if `TRUE`, enables tracing during optimization.
#' @return A list with optimized parameters.
optimize_this.gd <- function(setup, con, do_trace) {
  parameters <- con$initial_parameters
  objective <- setup$H
  grad <- setup$H_grad
  gamma <- con$gamma
  epsilon <- con$stopping_tolerance

  small_relative_ascent <- function(parameters, old_parameters) {
    objective_diff <- objective(parameters) - objective(old_parameters)
    objective_diff >= 0 &&
      objective_diff <= epsilon * (abs(objective(parameters)) + epsilon)
  }

  converged <- FALSE
  parameters <- parameters
  while (!converged) {
    old_parameters <- parameters
    gradient <- grad(old_parameters)
    parameters <- old_parameters + gamma * gradient
    if (small_relative_ascent(parameters, old_parameters)) {
      converged <- TRUE
    }
  }

  list(
    parameters = parameters,
    con = con,
    setup
  )
}

#' Maximum Likelihood Estimation
#'
#' Calculates the maximum likelihood estimators for `mu` and `sigma2`.
#'
#' @param x A numeric vector of data.
#' @param w A numeric vector of weights.
#' @param nu Degrees of freedom parameter.
#' @return A numeric vector with estimated `mu` and `sigma2`.
#' @export
calculate_mle <- function(x, w, nu) {
  n <- length(x)

  mu_hat <- sum(x * w) / sum(w)
  sigma2_hat <- sum(w * (x - mu_hat)^2) / (n * nu)

  c(mu = mu_hat, sigma2 = sigma2_hat)
}
