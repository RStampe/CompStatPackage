get_simulated_data <- function(n = 2000, beta = c(1, 0.5, -1, 0.5, 1)) {
  X <- runif(n, min = 0, max = 10)  # Uniform distribution between 0 and 10
  knots_inner <- get_inner_knots(X, length(beta))
  knots <- get_knots(knots_inner)
  Omega = calculate_Omega(knots_inner)
  X_spline <- splineDesign(X, knots = knots)
  num_basis <- ncol(X_spline)
  eta <- X_spline %*% beta  # No intercept
  p <- 1 / (1 + exp(-eta))
  Y <- rbinom(n, size = 1, prob = p)
  simulated_data <- list(x = list(x = X_spline, Omega = Omega), y = Y, beta = beta)
  print(summary(simulated_data))
  simulated_data
}


get_inner_knots <- function(X, length_beta) {
  seq(min(X), max(X), length.out = length_beta - 2)
}

get_knots <- function(inner_knots) {
  sort(c(rep(range(inner_knots), 3), inner_knots))
}

calculate_phi <- function(X, knots) {
  splineDesign(knots, X)
}

calculate_Omega <- function(inner_knots) {
  knots <- get_knots(inner_knots)
  d <- diff(inner_knots) # The vector of knot differences; b - a
  g_ab <- splineDesign(knots, inner_knots, derivs = 2)
  knots_mid <- inner_knots[-length(inner_knots)] + d / 2
  g_ab_mid <- splineDesign(knots, knots_mid, derivs = 2)
  g_a <- g_ab[-nrow(g_ab), ]
  g_b <- g_ab[-1, ]
  (crossprod(d * g_a, g_a) +
      4 * crossprod(d * g_ab_mid, g_ab_mid) +
      crossprod(d * g_b, g_b)) / 6
}

calculate_H <- function(X, y, Omega, lambda, beta) {
  f_x <- X %*% beta

  penalizer <- lambda * crossprod(beta, crossprod(Omega, beta))

  -mean(y * f_x - log(1 + exp(f_x))) + penalizer
}

calculate_grad_H_in_all_point <- function(X, y, Omega, lambda, beta) {
  exp_X_beta <- exp(X %*% beta)

  -t(crossprod((y - exp_X_beta / (1 + exp_X_beta)), X) / length(y)) + 2 * lambda * crossprod(Omega, beta)
}

calculate_grad_H <- function(X, y, Omega, lambda, beta, i) {
  n <- length(i)
  X <- X[i, , drop = FALSE]
  y <- y[i]
  Xbeta <- X %*% beta
  temp <- plogis(Xbeta)

  residuals <- y - temp
  gradient <- -1 / n * crossprod(X, residuals)

  if(lambda != 0) {
    gradient <- gradient + 2 * lambda * Omega %*% beta
  }

  gradient
}



calculate_Hessian_H_in_all_point <- function(X, y, Omega, lambda, beta) {

  s <- 1 / (1 + exp(-X %*% beta))  # N x 1 vector

  weight_vector <- as.vector(s * (1 - s))  # N x 1 vector

  # Compute X^T * W * X efficiently without forming the full W matrix
  X_weighted <- X * weight_vector  # Element-wise multiplication, N x p matrix
  Hessian <- (1 / length(y)) * t(X) %*% X_weighted + 2 * lambda * Omega  # p x p matrix

  return(Hessian)
}


get_design_do_nothing <-  function(x){
  if(!"Omega" %in% names(x)){
    warning("Omega is 0")
    Omega = diag(0,ncol(x))
  }
  else {
    Omega = x$Omega
  }

  if(!"x" %in% names(x)){
    design = x
  }
  else {
    design = x$x
  }


  structure(list(X = design, Omega = Omega), class = "design")
}

get_design_B_splines <- function(x_data, parameter_size) {
  inner_knots <- get_inner_knots(x_data, parameter_size)
  knots <- get_knots(inner_knots)
  X <- calculate_phi(x_data, knots)
  Omega <- calculate_Omega(inner_knots)

  structure(
    list(
      X = X,
      Omega = Omega,
      knots = knots
    ),
    class = "design"
  )
}


get_design_DR = function(x, parameter_size){
  initial <- get_design_B_splines(x, parameter_size)
  X = initial$X
  Omega = initial$Omega
  X_svd = svd(X)
  Omega_tilde = t(crossprod(X_svd$v, Omega %*% X_svd$v)) / X_svd$d
  Omega_tilde = t(Omega_tilde) / X_svd$d
  Omega_tilde_svd = svd(Omega_tilde)
  U_tilde = X_svd$u %*% Omega_tilde_svd$u

  structure(
    list(
      X = U_tilde,
      Omega = Omega_tilde,
      knots = initial$knots
    ),
    class = "design"
  )
}


function_H_factory <- function(design, y, calculate_grad_H = calculate_grad_H) {
  if (!inherits(design, "design")) {
    warning("design object is not a design class")
  }

  if (!all(c("X", "Omega") %in% names(design))) {
    stop("Missing one or both of X and Omega.")
  }
  X <- design$X
  Omega <- design$Omega

  get_h_objective <- function() {
    function(beta, lambda) calculate_H(X, y, Omega, lambda, beta)
  }

  get_gradient_H_all_points <- function() {
    function(beta, lambda) {
      calculate_grad_H_in_all_point(X, y, Omega, lambda, beta)
    }
  }

  get_gradient_H_i <- function() {
    function(beta, lambda, i) {
      calculate_grad_H(X, y, Omega, lambda, beta, i)
    }
  }

  get_Hessian <- function() {
    function(beta, lambda) {
      calculate_Hessian_H_in_all_point(X, y, Omega, lambda, beta)
    }
  }

  get_Hessian_numeric <- function() {
    function(beta, lambda) {
      hessian(function(beta) calculate_H(X, y, Omega, lambda, beta), beta)
    }
  }

  structure(list(
    objective = get_h_objective(),
    gradient = get_gradient_H_all_points(),
    gradient_i = get_gradient_H_i(),
    Hessian = get_Hessian(),
    Hessian_numeric = get_Hessian_numeric()
  ), class = "functions")
}


function_loglike_factory <- function(design, y) {


  if (!inherits(design, "design")) {
    warning("design object is not a design class")
  }

  if (!all(c("X") %in% names(design))) {
    stop("Missing one or both of X and Omega.")
  }
  X <- design$X

  # Negative log-likelihood function
  calculate_H <- function(x, y, beta) {
    residuals <- y - x %*% beta
    H <- 0.5 * sum(residuals^2)
    return(H)
  }

  # Gradient of the negative log-likelihood
  calculate_grad_H_in_all_point <- function(x, y, beta) {
    residuals <- y - x %*% beta
    grad_H <- -t(x) %*% residuals
    return(as.vector(grad_H))
  }

  # Gradient at a single data point
  calculate_grad_H <- function(x, y, beta, i) {
    x_i <- x[i, , drop = FALSE]
    y_i <- y[i]
    residual_i <- y_i - x_i %*% beta
    grad_H_i <- -t(x_i) %*% residual_i
    return(as.vector(grad_H_i))
  }

  # Hessian matrix of the negative log-likelihood
  calculate_Hessian_H_in_all_point <- function(x, y, beta) {
    Hessian_H <- t(x) %*% x
    return(Hessian_H)
  }

  # Objective function
  get_h_objective <- function() {
    function(beta, lambda) calculate_H(x, y, beta)
  }

  # Gradient function
  get_gradient_H_all_points <- function() {
    function(beta, lambda) {
      calculate_grad_H_in_all_point(x, y, beta)
    }
  }

  # Gradient function at a single data point
  get_gradient_H_i <- function() {
    function(beta, lambda, i) {
      calculate_grad_H(x, y, beta, i)
    }
  }

  # Hessian function
  get_Hessian <- function() {
    function(beta, lambda) {
      calculate_Hessian_H_in_all_point(x, y, beta)
    }
  }

  # Numerical Hessian using numDeriv package
  get_Hessian_numeric <- function() {
    function(beta, lambda) {
      numDeriv::hessian(function(b) calculate_H(x, y, b), beta)
    }
  }

  # Closed-form solution for MLE
  get_MLE_solution <- function() {
    beta_hat <- solve(t(x) %*% x, t(x) %*% y)
    return(as.vector(beta_hat))
  }

  structure(list(
    objective = get_h_objective(),
    gradient = get_gradient_H_all_points(),
    gradient_i = get_gradient_H_i(),
    Hessian = get_Hessian(),
    Hessian_numeric = get_Hessian_numeric(),
    MLE_solution = get_MLE_solution()
  ), class = "functions")
}

setup_for_optim <- function(X, y, function_factory, design_function, num_parameters = NULL, ...) {
  design_matrix <- if (is.null(num_parameters)) {
    design_function(X)
  } else {
    design_function(X, num_parameters)
  }

  model_functions <- function_factory(design_matrix, y, ...)

  if (!inherits(model_functions, "functions")) {
    stop("model_functions object is not a model_functions class")
  }

  structure(list(design = design_matrix, y = y, model_functions = model_functions), class = "setup")
}

perform_optim <- function(setup, algorithm, configuration, ...) {
  UseMethod("perform_optim")
}


#' Perform Optimization Based on Setup Configuration
#'
#' `perform_optim.setup` runs an optimization routine using the provided setup
#' @param setup list
#' @param algorithm algorithms.
#' @param configuration A list of configuration
#' @param ... Additional arguments
#' @return A list or object
#' @export
perform_optim.setup <- function(setup, algorithm, configuration, ...) {
  default_configuration <- list(
    lambda = 0,
    initial_parameters = rep(0, ncol(setup$design$X)),
    num_samples = length(setup$y),  # Dynamic based on setup
    learning_rate = decay_scheduler(),
    max_num_epochs = 250,
    mini_batch_size = 50,
    sampling_function = sample,
    batching = batch,  # For batching algorithm either batch, momentum or adam
    clipping = T,
    epsilon = 1e-5,
    patience = 5,
    stopping_criteria = stopping_objective
  )
  # Merge default configuration with user-provided configuration
  configuration <- modifyList(default_configuration, configuration)

  if(exists("use_stopping")) if(!use_stopping) configuration$stopping_criteria = function(...) F

  if(exists("my_debug")) if(my_debug) browser()

  algorithm(setup, configuration, ...)
}

#### Naive implementation ####
sgd_naive <- function(setup, configuration) {
  lambda <- configuration$lambda
  initial_parameters <- configuration$initial_parameters
  num_samples <- configuration$num_samples
  learning_rate <- configuration$learning_rate
  max_num_epochs <- configuration$max_num_epochs
  sampling_function <- configuration$sampling_function

  gradient_function <- setup$model_functions$gradient_i
  learning_rate <- if (is.function(learning_rate)) {
    learning_rate(current_epoch = 1:max_num_epochs)
  } else {
    rep(learning_rate, max_num_epochs)
  }

  for (current_epoch in 1:max_num_epochs) {
    sampled_indices_list <- sampling_function(num_samples)
    for (j in 1:num_samples) {
      index <- sampled_indices_list[j]
      initial_parameters <- initial_parameters - learning_rate[current_epoch] * gradient_function(initial_parameters, lambda, index)
    }
  }
  structure(
    list(
      optimized_parameters = initial_parameters,
      objective_value = setup$model_functions$objective(initial_parameters, lambda),
      lambda = lambda,
      model = setup
    ),
    class = "optim_result"
  )
}
#### Advanced implementation ####
decay_scheduler <- function(gamma0 = 1, a = 1, K = 1, gamma1, n1) {
  force(a)
  if (!missing(gamma1) && !missing(n1))
    K <- n1^a * gamma1 / (gamma0 - gamma1)
  b <- gamma0 * K
  function(n) b / (K + n^a)
}

clipped_gradient <- function(grad, clip_value = 1000) {
  ifelse((abs(grad) > clip_value), sign(grad) * clip_value, grad)
}


stopping_objective <- function(n, obj_function, current_parameter, con, best_loss_ref, no_improve_ref) {
  if (n %% 5 == 0) {  # Check every 5 epochs
    val_loss <- obj_function(current_parameter, con$lambda)

    if (val_loss < best_loss_ref - con$epsilon) {
      best_loss_ref <- val_loss  # Update best_loss in the calling scope
      no_improve_ref <- 0       # Reset no_improve count
    } else {
      no_improve_ref <- no_improve_ref + 1
    }

    # Early stopping based on patience
    if (no_improve_ref >= con$patience) {
      message("Early stopping at epoch ", n, " with best validation loss: ", best_loss_ref)
      return(c(1, best_loss_ref, no_improve_ref ))
    }
  }
  return(c(0, best_loss_ref, no_improve_ref ))
}


batch <- function(
    parameter,
    sample_indices,
    learning_rate,
    gradient_function,
    lambda,
    mini_batch_size,
    clipping,
    ...) {
  M <- floor(length(sample_indices) / mini_batch_size)
  for (j in 0:(M - 1)) {
    i <- sample_indices[(j * mini_batch_size + 1):(j * mini_batch_size + mini_batch_size)]
    parameter <- parameter - learning_rate * clipping(gradient_function(parameter, lambda, i))
  }
  parameter
}

adam <- function() {
  rho <- v <- 0
  function(
    parameter,
    sample_indices,
    learning_rate,
    gradient_function,
    lambda,
    mini_batch_size,
    adam_momentum_mermory_1 = 0.9,     # Momentum memory
    adam_momentum_mermory_2 = 0.9,     # Second moment memory
    ...
  ) {
    M <- floor(length(sample_indices) / mini_batch_size)
    for(j in 0:(M - 1)) {
      i <- sample_indices[(j * mini_batch_size + 1):(j * mini_batch_size + mini_batch_size)]
      gr <- gradient_function(parameter, lambda, i)
      rho <<- adam_momentum_mermory_1 * rho + (1 - adam_momentum_mermory_1) * gr
      v <<- adam_momentum_mermory_2 * v + (1 - adam_momentum_mermory_2) * gr^2
      parameter <- parameter - learning_rate * (rho / (sqrt(v) + 1e-8))
    }
    parameter
  }
}

momentum <- function() {
  rho <- 0
  function(
    parameter,
    sample_indices,
    learning_rate,
    gradient_function,
    lambda,
    mini_batch_size,
    momentum_mermory = 0.95,        # Momentum memory
    ...
  ) {
    M <- floor(length(sample_indices) / mini_batch_size)
    for(j in 0:(M - 1)) {
      i <- sample_indices[(j * mini_batch_size + 1):(j * mini_batch_size + mini_batch_size)]
      # Using '<<-' assigns the value to rho in the enclosing environment
      rho <<- momentum_mermory * rho + (1 - momentum_mermory) * gradient_function(parameter, lambda, i)
      parameter <- parameter - learning_rate * rho
    }
    parameter
  }
}

sgd_wrapper <- function(setup, con, ...) {

  # needed for stopping criteria
  best_loss <- Inf
  no_improve <- 0
  patience <- con$patience  # Number of epochs without improvement to stop
  stopping_function <- con$stopping_criteria

  clipping_function = if(con$clipping){clipped_gradient} else {function(grad) grad}

  current_parameter <- con$initial_parameters
  gradient_function <- setup$model_functions$gradient_i

  if (is.function(con$learning_rate)){
    learning_rate <- con$learning_rate(1:con$max_num_epochs)
  } else {
    learning_rate <- rep_len(con$learning_rate, con$max_num_epochs)
  }
  for (n in 1:con$max_num_epochs) {

    sample_indices <- con$sampling_function(con$num_samples)
    current_parameter <- con$batching(parameter = current_parameter,
                                      sample_indices = sample_indices,
                                      learning_rate = learning_rate[n],
                                      gradient_function = gradient_function,
                                      lambda = con$lambda,
                                      mini_batch_size = con$mini_batch_size,
                                      clipping = clipping_function,
                                      ...)


    res_stopping <- stopping_function(n,
                                      obj_function = setup$model_functions$objective,
                                      current_parameter,
                                      con,
                                      best_loss_ref = best_loss,
                                      no_improve_ref = no_improve)
    if(length(res_stopping) > 1) {
      best_loss <- res_stopping[2]
      no_improve <- res_stopping[3]
    }
    if (res_stopping[1]) break

  }
  structure(
    list(
      optimized_parameters = current_parameter,
      objective_value = setup$model_functions$objective(current_parameter, con$lambda),
      model = setup,
      learning_rate = learning_rate,
      con = con
    ),
    class = "optim_result"
  )
}


#### Other algorithm types ####
gradient_descent <- function(setup, con) {
  current_parameter <- con$initial_parameters




  backtracking_factor <-0.8  # For gradient descent
  armijo_constant <- 0.1      # For gradient descent
  convergence_threshold <- 1e-4 # For gradient descent
  # Is set by line search
  learning_rate = 1

  best_loss <- Inf
  no_improve <- 0
  patience <- con$patience  # Number of epochs without improvement to stop
  stopping_function <- con$stopping_criteria


  objective_function <- setup$model_functions$objective
  gradient_function <- setup$model_functions$gradient
  for (n in 1:con$max_num_epochs) {

    value <- objective_function(current_parameter, con$lambda)
    observed_gradient <- gradient_function(current_parameter, con$lambda)
    gradient_norm <- sum(observed_gradient^2)
    # Convergence criterion based on gradient norm
    # if (gradient_norm <= convergence_threshold) break
    learning_rate <- learning_rate

    # Proposed descent step
    updated_parameters <- current_parameter - learning_rate * observed_gradient
    # Backtracking while descent is insufficient
    while (objective_function(updated_parameters, con$lambda) > value - armijo_constant * learning_rate * gradient_norm) {
      learning_rate <- backtracking_factor * learning_rate


      updated_parameters <- current_parameter - learning_rate * observed_gradient
    }
    current_parameter <- updated_parameters


    res_stopping <- stopping_function(n,
                                      obj_function = setup$model_functions$objective,
                                      current_parameter,
                                      con,
                                      best_loss_ref = best_loss,
                                      no_improve_ref = no_improve)
    if(length(res_stopping) > 1) {
      best_loss <- res_stopping[2]
      no_improve <- res_stopping[3]
    }

    if (res_stopping[1]) break


  }
  if (n == con$max_num_epochs) {
    warning("Maximum number, ", con$max_num_epochs, ", of iterations reached")
  }
  structure(
    list(
      optimized_parameters = current_parameter,
      objective_value = setup$model_functions$objective(current_parameter, con$lambda),
      model = setup,
      learning_rate = learning_rate
    ),
    class = "optim_result"
  )
}

newton <- function(setup, con) {
  objective_function  <- setup$model_functions$objective
  gradient_function <- setup$model_functions$gradient
  hessian_function  <- setup$model_functions$Hessian

  current_parameter <- con$initial_parameters


  best_loss <- Inf
  no_improve <- 0
  patience <- con$patience  # Number of epochs without improvement to stop

  max_num_epochs <- con$max_num_epochs
  backtracking_factor = 0.8
  armijo_constant = 0.1
  gamma0 = 1
  epsilon = con$epsilon

  for(n in 1:max_num_epochs) {
    value <- objective_function(current_parameter, con$lambda)

    # Update early stopping criteria based on validation loss
    if (value < best_loss - epsilon) {  # Check if loss improves
      best_loss <- value
      no_improve <- 0  # Reset improvement counter
    } else {
      no_improve <- no_improve + 1
    }

    # Stop if no improvement after the specified patience interval
    if (no_improve >= patience) {
      message("Early stopping at epoch ", n, " with best validation loss: ", best_loss)
      break
    }

    grad <- gradient_function(current_parameter, con$lambda)
    grad_norm_sq <- sum(grad^2)

    Hessian <- hessian_function(current_parameter, con$lambda)
    rho <- - drop(solve(Hessian, grad))
    gamma <- gamma0
    par1 <- current_parameter + gamma * rho
    h_prime <- t(grad) %*% rho
    while(objective_function(par1, con$lambda) > value +  armijo_constant * gamma * h_prime) {
      gamma <- backtracking_factor * gamma
      par1 <- current_parameter + gamma * rho
    }
    current_parameter <- par1
  }
  if(n == max_num_epochs)
    warning("Maximum number, ", max_num_epochs, ", of iterations reached")

  structure(
    list(
      optimized_parameters = current_parameter,
      objective_value = setup$model_functions$objective(current_parameter, con$lambda),
      model = setup
    ),
    class = "optim_result"
  )
}



#' Optimization Function for Model Fitting
#'
#' `my_optim` performs optimization on input data and target values using
#' specified basis functions and regularization. It provides an option to
#' select between batch and other algorithms and allows custom control settings.
#'
#' @param X A matrix or data frame of predictor variables.
#' @param y A numeric vector of target values.
#' @param basis A character string specifying the basis function to "splines", "DR", "non".
#' @param p A numeric value representing parameters for the basis function if applicable.
#' @param lambda A numeric value specifying the regularization parameter.
#' @param algorithm A character string indicating the optimization algorithm "batch", "adam", "momentum", "newton", "gd".
#' @param con A list of control settings for the optimization process, such as lambda = 0, initial_parameters, num_samples, learning_rate, max_num_epochs, mini_batch_size, sampling_function, clipping, epsilon, patience, stopping_criteria.
#' @param ... Additional arguments to be passed to the optimization algorithm.
#'
#' @return A list containing the optimized parameters, final loss value,
#'   and convergence status.
#' @export
my_optim <- function(X, y, basis = "non", p, lambda, algorithm = "batch", con = list(),...) {
  basis <- match.arg(basis, c("splines", "DR", "non"))
  design_function <- ifelse(basis == "splines", get_design_B_splines, ifelse(basis == "DR", get_design_DR, get_design_do_nothing))
  setup <- setup_for_optim(X, y, function_H_factory, design_function, num_parameters = p)

  algorithm <- match.arg(algorithm, c("batch", "adam", "momentum", "newton", "gd"))
  algorithm_wrapper<- ifelse(algorithm == "newton", newton, ifelse(algorithm == "gd", gradient_descent, sgd_wrapper))

  con$batching <- ifelse(basis == "adam", adam, ifelse(basis == "momentum", momentum, batch))
  con$lambda = lambda
  perform_optim(setup, algorithm_wrapper, con)
}


