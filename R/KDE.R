#' Title: Epanechnikov Kernel Function
#'
#' Description: Calculates the value of the Epanechnikov kernel for a given input, `x`.
#' The Epanechnikov kernel is a common kernel function used in density estimation.
#'
#' @param x A numeric value for which the kernel value is to be calculated.
#' @return A numeric value representing the kernel density, which will be zero for absolute values of `x` greater than 1.
#' @export
epanechnikov_kernel <- function(x) (abs(x) <= 1) * (3 / 4) * (1 - x^2)


#' Title: Gaussian Kernel Function
#'
#' Description: Calculates the value of the Gaussian kernel for a given input, `x`.
#' The Gaussian kernel is widely used in kernel density estimation and smoothing.
#'
#' @param x A numeric value for which the Gaussian kernel density is to be calculated.
#' @return A numeric value representing the Gaussian kernel density at the input `x`.
#' @export
gaussian_kernel <- function(x) (1 / sqrt(2 * pi)) * exp(-x^2 / 2)

#' Title: Kernel Density Estimator for a Vector of Grid Points
#'
#' Description: Computes kernel density estimates at specified grid points using a given kernel function.
#'
#' @param grid_points A numeric vector of points where density estimates are evaluated.
#' @param vec_obs A numeric vector of observed data points.
#' @param bandwidth A positive numeric value representing the smoothing bandwidth.
#' @param kernel_function A function specifying the kernel to use, which should take a numeric vector as input.
#' @return A list with two components: `x` (grid points) and `y` (corresponding density estimates).
#' @export
kernel_vec <- function(grid_points, vec_obs, bandwidth, kernel_function) {
  vec_density_values <- numeric(length(grid_points))

  for (i in seq_along(grid_points)) {
    x <- (grid_points[i] - vec_obs) / bandwidth
    vec_density_values[i] <- sum(kernel_function(x)) / (bandwidth * length(vec_obs))
  }

  list(x = grid_points, y = vec_density_values)
}


#' Title: Optimized Kernel Density Estimator with Binning
#'
#' Description: Calculates kernel density estimates at specified grid points using binning to improve efficiency.
#'
#' @param grid_points A numeric vector of grid points where density estimates are evaluated.
#' @param vec_obs A numeric vector of observed data points.
#' @param bandwidth A positive numeric value representing the smoothing bandwidth.
#' @param kernel_function A kernel function that takes a numeric vector as input.
#' @return A list with components: `x` (grid points) and `y` (corresponding density estimates).
#' @export
kernel_binning_smart <- function(grid_points, vec_obs, bandwidth, kernel_function) {
  grid_lower <- min(grid_points)
  grid_upper <- max(grid_points)
  grid_length <- length(grid_points)
  grid_diff <- grid_points[2] - grid_lower
  max_nonzero_i <- floor(bandwidth / grid_diff)
  kernel_evals <- kernel_function((grid_points[1:(max_nonzero_i + 1)] - grid_lower) / bandwidth)
  kernel_vec <- c(rev(kernel_evals[-1]), kernel_evals)
  weights <- c(rep(0, max_nonzero_i),
               binning_cpp(vec_obs, grid_lower, grid_upper, grid_length),
               rep(0, max_nonzero_i))
  y <- numeric(grid_length)
  for (i in (1 + max_nonzero_i):(grid_length + max_nonzero_i))
    y[i - max_nonzero_i] <- sum(weights[(i - max_nonzero_i):
                                          (i + max_nonzero_i)] *
                                  kernel_vec)

  list(x = grid_points, y = y/bandwidth)
}


#' Title: Kernel Density Estimation
#'
#' Description: Computes a kernel density estimate for a given set of observations using a specified kernel function.
#' Offers flexibility in choosing a custom kernel function, kernel calculator, and number of grid points for estimation.
#'
#' @param vec_obs A numeric vector of observed data points.
#' @param bandwidth A positive numeric value representing the smoothing bandwidth.
#' @param kernel_function A kernel function (e.g., Gaussian, Epanechnikov) that takes a numeric vector as input.
#'                        Defaults to Gaussian if not specified.
#' @param kernel_calculator Optional alternative kernel density calculator function for custom density calculations.
#' @param n The number of grid points to use for density estimation (default is 512).
#' @param ... Additional arguments to pass to `kernel_function` or `kernel_calculator`.
#' @return A list containing:
#'   - `x`: Grid points where the density is estimated.
#'   - `y`: Corresponding density values.
#' @export
kernel_density_estimation <- function(vec_obs, bandwidth, kernel_function = NULL, kernel_calculator = NULL, n = 512, ...) {

  grid_points <- seq(min(vec_obs) - 3 * bandwidth, max(vec_obs) + 3 * bandwidth, length.out = n)

  if (is.null(kernel_function)) {
    if (grepl("epanechnikov", deparse(substitute(kernel_calculator)))) bandwidth <- bandwidth * sqrt(5)

    xy_values <- kernel_calculator(grid_points, vec_obs, bandwidth)
  } else {
    if (grepl("epanechnikov", deparse(substitute(kernel_function)))) bandwidth <- bandwidth * sqrt(5)

    xy_values <- kernel_calculator(grid_points, vec_obs, bandwidth, kernel_function, ...)
  }

  result <- list(
    x = xy_values$x, # Grid points (x-axis)
    y = xy_values$y, # Density values (y-axis)
    bw = bandwidth, # Bandwidth used
    n = length(vec_obs), # Number of observations
    call = match.call(), # Call to the function (similar to density())
    data = vec_obs,
    data.name = deparse(substitute(vec_obs)), # Name of the data
    has.na = anyNA(vec_obs)
  ) # Flag for missing values

  # Set the class to "density"
  class(result) <- c("density", "my_density")

  return(result)
}

#' Title: Silverman's Rule of Thumb for Bandwidth
#'
#' Description: Computes the bandwidth for kernel density estimation using Silverman's rule of thumb.
#'
#' @param vec_obs A numeric vector of observed data points.
#' @return A numeric value representing the estimated bandwidth.
#' @export
bandwidth_silvermans_rule <- function(vec_obs) 0.9 * min(var(vec_obs), IQR(vec_obs) / 1.34) * length(vec_obs)^(-1 / 5)

#' Title: Oracle Bandwidth for Epanechnikov Kernel
#'
#' Description: Computes an optimal (oracle) bandwidth for kernel density estimation with the Epanechnikov kernel.
#'
#' @param vec_obs A numeric vector of observed data points.
#' @return A numeric value representing the oracle bandwidth for the Epanechnikov kernel.
#' @export
calculate_oracle_bw_epanechnikov <- function(vec_obs) {

  evaluate_double_sum_of_kernel_vec <- function(vec_obs, pilot_bw) {
    n <- length(vec_obs)

    d4_std_gaussian_density <- function(x_i, x_j, pilot_bw) {
      x <- (x_i - x_j) / pilot_bw
      exp(-x^2 / 2) * ((x^2 - 6) * x^2 + 3) / sqrt(2 * pi)
    }

    result <- 0
    for (i in 1:n) {
      result <- result + sum(d4_std_gaussian_density(vec_obs[i], vec_obs, pilot_bw))
    }

    result
  }


  pilot_bandwidth <- bandwidth_silvermans_rule(vec_obs)

  pilot_bandwidth <- sqrt(2) * pilot_bandwidth

  n <- length(vec_obs)

  double_sum <- evaluate_double_sum_of_kernel_vec(vec_obs, pilot_bandwidth)

  d2_f_norm_sq <- 1 / (n^2 * (pilot_bandwidth)^5) * double_sum


  # the variance sigma_K^2 of a epanechnikov distribution is 1/5
  sigma_K <- 0.4472135955
  sigma2_K <- 0.2
  sigma4_K <- 0.4
  # the norm of a epanechnikov kernel, int_-1^1 (K(x))^2 dx = 0.6
  K_norm_sq <- 0.6


  n <- length(vec_obs)

 (K_norm_sq / (n * sigma4_K * d2_f_norm_sq))^0.2
}


#' Title: Cross-Validation Bandwidth Selection
#'
#' Description: Calculates the optimal bandwidth for kernel density estimation using cross-validation.
#'
#' @param vec_obs A numeric vector of observed data points.
#' @param number_of_partitions The number of partitions for cross-validation (default is 10).
#' @param kernel_evaluator A kernel function, typically an optimized C++ version for efficiency.
#' @return A numeric value representing the CV-selected bandwidth.
#' @export
get_cv_bw <- function(vec_obs, number_of_partitions = 10, kernel_evaluator = kernel_cpp_epanechnikov) {

  get_cv_objective <- function(vec_obs, number_of_partitions, kernel_evaluator) {
    n <- length(vec_obs)

    loglike <- 0

    indicies <- rep(1:number_of_partitions, times = ceiling(n / number_of_partitions))[sample(n)]

    objective_function <- function(bandwidth) {
      for (i in 1:number_of_partitions) {
        x_j <- vec_obs[indicies != i]
        x_i <- vec_obs[indicies == i]

        list_of_point <- kernel_evaluator(
          grid_points = x_i,
          vec_obs = x_j,
          bandwidth = bandwidth
        )

        loglike <- loglike + sum(log(list_of_point$y))
      }
      if (is.infinite(loglike)) {
        loglike <- -1e8
      }
      loglike
    }

    objective_function
  }

  objective <- get_cv_objective(vec_obs, number_of_partitions, kernel_evaluator)

  optimize(
    f = objective,
    interval = c(0, max(vec_obs) - min(vec_obs)),
    # interval = c(0, 1.144 * sqrt(var(vec_obs)) * length(vec_obs)^(-1/5)),
    maximum = TRUE
  )$maximum
}


#' Title: Cross-Validation Bandwidth Selection
#'
#' Description: Calculates the optimal bandwidth for kernel density estimation using cross-validation.
#'
#' @param vec_obs A numeric vector of observed data points.
#' @param kernel_evaluator A kernel function, typically an optimized C++ version for efficiency.
#' @return A numeric value representing the CV-selected bandwidth.
#' @export
get_ucv_bw <- function(vec_obs, kernel_evaluator = kernel_cpp_epanechnikov) {

  get_UCV_objective <- function(vec_obs, kernel_evaluator) {
    if (typeof(vec_obs) != "double") {
      stop("vec_obs argument is not numeric")
    }
    n <- length(vec_obs)

    max_obs <- max(vec_obs)
    min_obs <- min(vec_obs)

    # we perform leave one out cv therefore number of partitions is set to n
    number_of_partitions <- n

    indicies <- rep(1:number_of_partitions, times = ceiling(n / number_of_partitions))[sample(n)]

    objective_function <- function(bandwidth) {
      est_mean_f_hat <- 0

      for (i in 1:number_of_partitions) {
        x_j <- vec_obs[indicies != i]

        x_i <- vec_obs[indicies == i]

        val <- sum(kernel_evaluator(
          grid_points = x_i,
          vec_obs = x_j,
          bandwidth = bandwidth)$y)


        est_mean_f_hat <- val + est_mean_f_hat
      }

      f_hat_sq <- function(x) {
        (kernel_evaluator(x, vec_obs, bandwidth)$y)^2
      }

      norm <- integrate(f_hat_sq, lower = min_obs, upper = max_obs)$value

      norm - 2 * est_mean_f_hat / number_of_partitions
    }
    objective_function
  }

  objective <- get_UCV_objective(vec_obs, kernel_evaluator)

  optimize(
    f = objective,
    interval = c(0, max(vec_obs) - min(vec_obs)),
    # interval = c(0, 1.144 * sqrt(var(vec_obs)) * length(vec_obs)^(-1/5)),
    maximum = FALSE
  )$minimum
}


#' Title: Kernel Density Estimation with Customizable Bandwidth and Kernel
#'
#' Description: Computes a kernel density estimate for a given data vector with options for bandwidth selection,
#' kernel type, binning, and grid point resolution.
#'
#' @param vec_data A numeric vector of observed data points.
#' @param bandwidth A bandwidth value (numeric) or method (character: "silverman" or "cv" for cross-validation).
#'                  Defaults to "silverman".
#' @param kernel The kernel function to use for density estimation. Options are "epanechnikov" (default) or "gaussian".
#' @param binning Logical; if TRUE, uses binning to improve computational efficiency.
#' @param number_of_gridpoints The number of grid points used for density estimation (default is 512).
#' @return A list with:
#'   - `x`: Grid points where the density is estimated.
#'   - `y`: Corresponding density values.
my_density <- function(vec_data, bandwidth = "silverman", kernel = "epanechnikov", binning = TRUE, number_of_gridpoints = 512) {
  number_of_gridpoints <- 2^round(log(number_of_gridpoints, base = 2))

  # Process the kernel parameter
  if (is.character(kernel)) {
    kernel <- match.arg(kernel, choices = c("gaussian", "epanechnikov"))
    kernel_calculator <- switch(kernel,
                                "gaussian" = if(binning) {kernel_cpp_gaussian_binning} else {kernel_cpp_gaussian},
                                "epanechnikov" = if(binning) {kernel_cpp_epanechnikov_binning} else {kernel_cpp_epanechnikov},
                                stop("Unknown kernel"))
    kernel_function <- NULL
  } else if (is.function(kernel)) {
    kernel_calculator <- if(binning) {kernel_binning_smart} else {kernel_vec}
    kernel_function <- kernel
  } else {
    stop("Invalid kernel parameter")
  }

  # Process the bandwidth parameter
  if (is.character(bandwidth)) {
    bandwidth <- match.arg(bandwidth, choices = c("silverman", "oracle", "cv", "loocv", "ucv"))
    kernel_calculator_bw = if(is.character(kernel)){
      if(kernel == "epanechnikov"){
        kernel_cpp_epanechnikov
      }
    } else
    {
      kernel_cpp_epanechnikov
    }

    bw <- switch(bandwidth,
                 "silverman" = bandwidth_silvermans_rule(vec_data),
                 "oracle" = calculate_oracle_bw_epanechnikov(vec_data),
                 "cv" = get_cv_bw(vec_data, number_of_partitions = 10, kernel_evaluator = kernel_calculator_bw),
                 "loocv" = get_cv_bw(vec_data, number_of_partitions = length(vec_data), kernel_evaluator = kernel_calculator_bw),
                 "ucv" = get_ucv_bw(vec_data, kernel_evaluator = kernel_calculator_bw),
                 stop("Unknown bandwidth method")
    )
  } else if (is.function(bandwidth)) {
    bw <- bandwidth(vec_data)
  } else if (is.numeric(bandwidth)) {
    bw <- bandwidth
  } else {
    stop("Invalid bandwidth parameter")
  }

  # Perform kernel density estimation
  kernel_density_estimation(
    vec_obs = vec_data,
    bandwidth = bw,
    kernel_function = NULL,
    kernel_calculator = kernel_calculator,
    n = number_of_gridpoints
  )
}

