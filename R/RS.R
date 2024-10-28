#' Import Poisson Data and Set Data Vectors
#'
#' This function reads a CSV file containing Poisson data, extracts relevant vectors,
#' and sets them as global variables in C++ for further processing.
#'
#' @return None. Sets global variables `x` and `z` in C++.
#' @export
import_poisson_data <- function() {
  # Locate the CSV file within the package
  poisson_data_path <- system.file("extdata", "poisson.csv", package = "CompStatPackage")

  # Read the data
  poisson_data <- read.csv(file = poisson_data_path)

  # Set `x` and `z` in the global environment
  x <<- poisson_data$x
  z <<- poisson_data$z

  # Call the C++ function to set global vectors
  set_data_vectors(poisson_data$x, poisson_data$z)
  invisible(NULL)
}

#' Logarithm of Target Density Function
#'
#' Computes the logarithm of the target density function for Poisson-distributed data,
#' with vectorized operations over an input vector \code{y}.
#'
#' @param y A numeric vector input for evaluating the target density.
#' @return The log target density for each value of \code{y}.
#' @export
log_target_density_vectorize <- Vectorize(function(y) {
  ifelse((y>=0), sum(y*z*x-exp(y*x)), -Inf)
})


#' Derivative of Log Target Density Function
#'
#' Computes the derivative of the log target density function for a Poisson model.
#'
#' @param y A numeric input vector for derivative evaluation.
#' @return The derivative of the log target density for each \code{y}.
#' @export
d_log_target_density_vectorize <- Vectorize(function(y) {
  sum(z*x - (x * exp(y * x)))
})


#' Log Unnormalized Gaussian Density
#'
#' This function calculates the unnormalized logarithm of a Gaussian density given
#' a mean and standard deviation.
#'
#' @param x Numeric value to evaluate the density at.
#' @param mean The mean of the Gaussian distribution.
#' @param sd The standard deviation of the Gaussian distribution.
#' @return The log of the unnormalized Gaussian density.
#' @export
log_unnormalized_gaussian_density <- function(x, mean, sd) {
  -(x-mean)^2/(2*sd^2)
}


get_log_alpha_prime <- function(log_target_density, log_envelope_density) {

  diff_log_density <- function(x) log_envelope_density(x) - log_target_density(x)

  log_alpha_prime <- optimize(f = diff_log_density, interval = c(0, 1), maximum = FALSE)

  min(log_alpha_prime$objective, diff_log_density(0))
}


accept_prob_est_gauss = function(mu, sd, log_target_density){
  grid <- seq(0,1,length.out = 50)

  log_p <- function(x) log_unnormalized_gaussian_density(x, mu, sd)
  log_alpha_prime = get_log_alpha_prime(log_target_density, log_p)
  #normalization_constant_target <- 1.050152e-41

  # 1/sqrt(2*pi) = 0.3989423
  #normalization_constant_density <- 0.3989423/sd
  normalization <- 4.189501e-42/sd
  alpha <- exp(log_alpha_prime) * normalization
  as.numeric(!any(log_target_density(grid) + log_alpha_prime >= log_p(grid))) * alpha
}


get_optimal_parameters <- function(log_target_density) {

  objective_function <- function(parameters){
    mu <- parameters[1]
    sd <- parameters[2]
    -accept_prob_est_gauss(mu, sd, log_target_density_cpp)
  }
  # Use optim() to find the minimum of the function
  result <- optim(par = c(0, 1),  # Initial values for x and y
                  fn = objective_function,  # The objective function
                  method = "BFGS")  # Optimization method
  if(!is.null(result$message))
    warning(result$message)

  result
}


#' Gaussian Envelope Validation
#'
#' Checks if a proposed envelope function is valid, ensuring it bounds the target density.
#'
#' @param log_target_density A function to compute the log target density.
#' @param log_envelope_density A function to compute the log envelope density.
#' @param lower Lower bound of the interval to check.
#' @param upper Upper bound of the interval to check.
#' @param n Number of points to sample within the interval.
#' @return \code{TRUE} if the envelope is valid; \code{FALSE} otherwise.
#' @export
is_an_evnelope <- function(log_target_density, log_envelope_density, lower = 0, upper = 10, n = 1000) {
  grid <- seq(lower,upper,length.out = n)

  !any(log_target_density(grid) > log_envelope_density(grid))
}


#' Construct Gaussian Envelope
#'
#' Builds a Gaussian envelope around a target density with customizable options.
#'
#' @param log_target_density Function for computing the log target density.
#' @param type Character string specifying envelope type, either \code{"tight"} or \code{"manual"}.
#' @param mu Optional mean for the Gaussian envelope.
#' @param sd Optional standard deviation for the Gaussian envelope.
#' @return An envelope object with sampling and density properties.
#' @export
get_gaussian_envelope <- function(log_target_density, type = "manual", mu = 0, sd = 1) {

  # Change these
  type <- match.arg(type, choices = c("tight", "manual"))

  if(type == "manual"){
    est_alpha <- accept_prob_est_gauss(mu, sd, log_target_density)
  } else if (type == "tight") {
    optim_object <- get_optimal_parameters(log_target_density)

    mu <- optim_object$par[1]
    sd <- optim_object$par[2]
    est_alpha <- -optim_object$value
  }

  log_envelope_density <- function(x) log_unnormalized_gaussian_density(x, mean = mu, sd = sd)

  if(!is_an_evnelope(log_target_density, log_envelope_density))
    warning("proposed evnelope is not working")

  # Calculate alpha_prime by maximizing the ratio of target_density/envelope_density
  log_alpha_prime <- get_log_alpha_prime(log_target_density, log_envelope_density)

  # Create the envelope object as a list
  setup <- list(
    sampler = function(n) rnorm(n, mean = mu, sd = sd),
    log_target_density = log_target_density,
    log_envelope_density = function(x) log_envelope_density(x)-log_alpha_prime, #envelope: log(g(x)) + log(alpha) for log(f(x))
    log_alpha_prime = log_alpha_prime, # Normalization constant log_alpha_prime
    mu = mu,
    sd = sd,
    est_alpha = est_alpha
  )

  class(setup) <- c("envelope_class", "gaussian_envelope_class")

  return(setup)
}

get_proposals <- function(x, n, ...) {
  UseMethod("get_proposals")
}

get_accepts <- function(x, proposals, ...) {
  UseMethod("get_accepts")
}

controller <- function(x, sample_size, ...) {
  UseMethod("controller")
}

add_class <- function(x, new_class) {
  class(x) <- c(class(x), new_class)
  x
}


#' Get proposals vectorize
#' @export
get_proposals.vectorize <- function(setup, n) {
  setup$sampler(n)
}

#' Get accepts vectorize
#' @export
get_accepts.vectorize <- function(setup, proposals) {
  (runif(length(proposals)) <= exp(setup$log_target_density(proposals) - setup$log_envelope_density(proposals)))
}



#' controller vectorize
#' @export
controller.vectorize <- function(setup, sample_size, scaling_factor = 1, min_proposals = 1) {

  samples <- numeric(sample_size)  # Preallocate for efficiency
  accepted_samples <- 0            # Counter for accepted samples
  iteration <- 0                   # Track number of iterations
  sample_list <- list()             # Store all proposed samples
  total_attempts <- 0              # Track total attempts

  while (accepted_samples < sample_size) {
    iteration <- iteration + 1

    # Adjust number of proposals dynamically
    num_proposals <- floor(max(scaling_factor * (sample_size - accepted_samples), min_proposals))

    # Generate a batch of proposals using the sampler from the setup
    proposals <- get_proposals(setup, num_proposals)
    accepts <- get_accepts(setup, proposals)

    # Append accepted proposals to the sample list
    sample_list[[iteration]] <- proposals[accepts]
    accepted_samples <- accepted_samples + sum(accepts)

    # Update scaling factor based on acceptance rate from the first iteration
    total_attempts <- total_attempts + num_proposals
    if (iteration == 1) {
      scaling_factor <- scaling_factor * sample_size / sum(accepts)
    }
  }

  # Return the first 'sample_size' samples, flattened into a vector
  final_samples <- unlist(sample_list)[1:sample_size]

  # Compute acceptance rate
  acceptance_rate <- accepted_samples / total_attempts

  result <- list(
    samples = final_samples,
    acceptance_rate = acceptance_rate,
    log_target_density = setup$log_target_density
  )

  class(result) <- "rejection_sampler"

  result
}

#' get proposals cpp
#' @export
get_proposals.cpp <- function(x, n) {
  get_proposals.vectorize(x, n)
}


#' get accepts cpp
#' @export
get_accepts.cpp <- function(setup, proposals) {
  get_accepts_gaussian_cpp(setup$mu,
                           setup$sd,
                           setup$log_alpha_prime,
                           proposals,
                           runif(length(proposals)))
}



#' controller cpp
#' @export
controller.cpp <- function(setup, sample_size, scaling_factor = 1, min_proposals = 1) {
  if (!("gaussian_envelope_class" %in% class(setup))) {
    stop("Provide an gaussian envelope class.")
  }
  controller.vectorize(setup, sample_size, scaling_factor, min_proposals)
}


rejection_sampler_wrapper <- function(setup, sample_size, type = "vectorize") {
  if (!("envelope_class" %in% class(setup))) {
    stop("Provide an envelope class.")
  }

  setup <- add_class(setup, match.arg(type, c("cpp", "vectorize")))
  controller(setup, sample_size)
}


create_sampler_function <- function(n, c, Q, a, b, z) {
  cu <- c * runif(n)
  i <- findInterval(cu, Q)
  log(a[i] * exp(-b[i]) * (cu - Q[i]) + exp(a[i] * z[i])) / a[i]
}



create_log_envelope_density <- function(x, z, a, b){
  i <- findInterval(x, z)
  a[i] * x + b[i]
}



#' Create a Log Affine Envelope Function
#'
#' This function constructs a log affine envelope for a target density, based on a specified range (`start`, `end`), knots,
#' and optional custom functions for creating a sampler and computing the log envelope density.
#' The envelope can be used for importance sampling or other density estimation methods.
#'
#' @param log_target_density A function that computes the log of the target density at a given point.
#' @param d_log_target_density An optional function that computes the derivative of the log target density. Default is \code{NULL}.
#' @param knots A numeric vector of knot points that define the intervals for the log affine envelope. Default is \code{NULL}.
#' @param create_sampler_function A function for generating samples from the envelope. Defaults to \code{create_sampler_function_cpp}.
#' @param create_log_envelope_density A function for computing the log envelope density for given values. Defaults to \code{create_log_envelope_density_cpp}.
#' @param start A numeric value indicating the start of the range for the envelope. Default is \code{0}.
#' @param end A numeric value indicating the end of the range for the envelope. Default is \code{1}.
#' @return A list containing functions and parameters for the log affine envelope, which can be used to sample and evaluate densities.
#' @export
get_log_affine_envelope <- function(log_target_density,
                                    d_log_target_density = NULL,
                                    knots = NULL,
                                    create_sampler_function = create_sampler_function_cpp,
                                    create_log_envelope_density = create_log_envelope_density_cpp,
                                    start = 0,
                                    end = 1) {

  if (is.null(d_log_target_density))
    a <- numDeriv::grad(log_target_density, knots, method = "simple")
  else
    a <- d_log_target_density(knots)

  b <- log_target_density(knots) - a * knots


  n <- length(a)
  z <- c(start, (b[2:n] - b[1:(n - 1)]) / (a[1:(n - 1)] - a[2:n]), end)
  Q_cumsum <- cumsum(c(0, exp(b[1:n]) * (exp(a[1:n] * z[2:(n + 1)]) - exp(a[1:n] * z[1:n])) / a[1:n]))

  Q = Q_cumsum[1:n]
  c = Q_cumsum[n + 1]


  # Create functions need
  sampler_function <- function(n) {
    create_sampler_function(n, c, Q, a, b, z)
  }

  log_envelope_density <- function(x) {
    create_log_envelope_density(x, z, a, b)
  }

  # Create the envelope object as a list
  setup <- list(
    sampler = sampler_function,
    log_target_density = log_target_density,
    log_envelope_density = log_envelope_density, #envelope: log(g(x)) + log(alpha) for log(f(x))
    log_alpha_prime = 0, # Normalization constant log_alpha_prime
    est_alpha = NA
  )

  if(!is_an_evnelope(log_target_density, log_envelope_density, lower = start, upper = end-0.1))
    warning("proposed evnelope might not work")



  class(setup) <- c("envelope_class", "log_affine_envelope_class")

  setup
}



#' Perform Adaptive One-Time Placement
#'
#' Adjusts knot placement in the envelope function based on initial samples from rejection sampling.
#'
#' @param initial_envelope Initial envelope object.
#' @param log_target_density Function for computing the log target density.
#' @param d_log_target_density Function for computing the derivative of the log target density.
#' @param initial_samples Number of initial samples.
#' @param num_knots Number of knots to place adaptively.
#' @return A list containing the modified envelope, knot positions, and initial samples.
#' @export
adaptive_one_time_placement <- function(initial_envelope, log_target_density, d_log_target_density, initial_samples = 1000, num_knots = 5, ...) {

  initial_res <- rejection_sampler_wrapper(initial_envelope, initial_samples, ...)


  probs <- seq(0, 1, length.out = num_knots + 2)[-c(1, num_knots + 2)]  # Exclude 0% and 100%
  knots <- quantile(initial_res$samples, probs = probs)

  new_envelope <- get_log_affine_envelope(log_target_density, d_log_target_density, knots)

  return(list(
    envelope = new_envelope,
    knots = knots,
    samples = initial_res$samples
  ))
}


#' Adaptive Knot Placement for Log Affine Envelope
#'
#' This function adapts the placement of knots in the log affine envelope of a target density function.
#' The adaptive placement improves sampling efficiency by adjusting knot positions based on gradient information.
#'
#' @param initial_knots A numeric vector of initial knot positions for constructing the log affine envelope.
#' @param log_target_density A function that computes the log of the target density at a given point.
#' @param d_log_target_density A function that computes the derivative of the log target density.
#' @param sample_size An integer specifying the number of samples to use in each iteration. Default is \code{1000}.
#' @param max_iterations An integer specifying the maximum number of iterations for the adaptive algorithm. Default is \code{10}.
#' @param num_tries An integer specifying the number of attempts to adjust each knot position. Default is \code{10}.
#' @return A numeric vector of adaptively placed knots for constructing the log affine envelope.
#' @export
adaptive_placement <- function(initial_knots, log_target_density, d_log_target_density, sample_size = 1000, max_iterations = 10, num_tries = 10) {
  # initial_knots: Numeric vector of initial knot positions
  # sample_size: Integer, number of samples to draw in each iteration
  # max_iterations: Integer, maximum number of iterations to perform

  current_knots <- initial_knots
  best_time <- Inf
  best_envelope <- NULL
  best_knots <- NULL
  acceptance_rates <- c()
  total_times <- c()
  samples <- c()
  tried_num_knots <- c()
  for (iteration in 1:max_iterations) {
    # Create envelope with current knots
    envelope <- get_log_affine_envelope(log_target_density, d_log_target_density, current_knots)

    # Measure time to sample
    times <- c()
    for(i in 1:num_tries)
    {
      start_time <- Sys.time()
      res <- rejection_sampler_wrapper(envelope, sample_size)
      end_time <- Sys.time()
      times <- c(times, as.numeric(difftime(end_time, start_time, units = "secs")))
    }
    total_time <- median(times)


    acceptance_rates <- c(acceptance_rates, res$acceptance_rate)
    total_times <- c(total_times, total_time)
    tried_num_knots <- c(tried_num_knots, length(current_knots))

    if (total_time < best_time) {
      best_time <- total_time
      best_envelope <- envelope
      best_knots <- current_knots
    }

    # Adjust the number of knots for the next iteration
    num_knots <- length(current_knots) + 1  # Increment number of knots
    probs <- seq(0, 1, length.out = num_knots + 2)[-c(1, num_knots + 2)]  # Exclude 0% and 100%
    current_knots <- quantile(res$samples, probs = probs)
  }

  return(list(
    envelope = best_envelope,
    knots = best_knots,
    acceptance_rates = acceptance_rates,
    total_times = total_times,
    num_knots = tried_num_knots
  ))
}


#' Perform Rejection Sampling with Specified Envelope Type
#'
#' This function performs rejection sampling to generate samples from a target density, using either a Gaussian
#' or log affine envelope. The function allows customization of the envelope type and initial knots for adaptive
#' placement within the envelope.
#'
#' @param n An integer specifying the number of samples to generate.
#' @param log_target_density A function that computes the log of the target density at a given point.
#' @param envelope_type A character string specifying the type of envelope to use. Options are \code{"gaussian"} (default) or \code{"log_affine"}.
#' @param d_log_target_density An optional function that computes the derivative of the log target density. Default is \code{NULL}.
#' @param initial_knots A numeric vector of initial knot positions, used when \code{envelope_type} is \code{"log_affine"}. Default is \code{c(0.1, 0.5)}.
#' @return A numeric vector of accepted samples from the target density.
#' @export
perform_rejection_sampling <- function(n, log_target_density, envelope_type = "gaussian", d_log_target_density = NULL, initial_knots =  c(0.1,0.5)) {

  type <- match.arg(envelope_type, choices = c("gaussian", "log_affine"))

  if(type == "gaussian"){
    envelope <- get_gaussian_envelope(log_target_density, type = "tight")
  } else if (type == "log_affine") {

    if(n < 1e6) {
      envelope <- get_log_affine_envelope(log_target_density, d_log_target_density, knots = initial_knots)
      if(n >= 1e5)
        envelope <- adaptive_one_time_placement(envelope, log_target_density_cpp, d_log_target_density_cpp, num_knots = 2)$envelope
    } else {
      envelope <- adaptive_placement(
        initial_knots = initial_knots,
        log_target_density_cpp,
        d_log_target_density_cpp,
        sample_size = n/100,
        num_tries = 5,
        max_iterations = floor(log(n))
      )$envelope
    }
  }

  if ("gaussian_envelope_class" %in% class(envelope)) {
    type = "cpp"
  }
  else {
    type = "vectorize"
  }
  rejection_sampler_wrapper(envelope, sample_size = n, type)
}




plot_target <- function(object, ...) {
  UseMethod("plot_target")  # Tells R to dispatch methods based on the object's class
}



#' Plot Target Density
#'
#' Plots the target density of an envelope class object, optionally on a logarithmic scale.
#'
#' @param object An object of class \code{envelope_class} containing the target density function.
#' @param width A numeric value defining the additional range width for the plot beyond the peak point. Default is \code{1}.
#' @param log_scale A logical value indicating if the plot should use a logarithmic scale. Default is \code{TRUE}.
#' @param lower, upper Numeric values defining the lower and upper bounds for the optimization range. Default is \code{0} and \code{0.5}.
#' @return A ggplot object showing the target density plot.
#' @export
plot_target.envelope_class <- function(object, width = 1, log_scale = T, lower = 0, upper = 0.5) {

  top_point <- optimise(object$log_target_density, interval = c(lower, upper), maximum = TRUE)$maximum

  gridpoints <- seq(0, top_point+width, length.out = 512)

  transform <- function(x){
    if (!log_scale) {
      return(exp(x))
    }
    x
  }

  df <- data.frame(
    gridpoints = gridpoints,
    target = object$log_target_density(gridpoints)
  )

  # Plot the function
  ggplot(df, aes(x = gridpoints, y = transform(target))) +
    geom_line() +
    labs(title = "Plot of density", x = "y", y = "")
}

#' Plot Target and Envelope Densities
#'
#' Plots the target and envelope densities for an envelope class object, with options for logarithmic scale.
#'
#' @param object An object of class \code{envelope_class} containing the target and envelope density functions.
#' @param width A numeric value defining the additional range width for the plot beyond the peak point. Default is \code{1}.
#' @param log_scale A logical value indicating if the plot should use a logarithmic scale. Default is \code{TRUE}.
#' @param lower, upper Numeric values defining the bounds for optimization and density evaluation. Default is \code{0} and \code{0.5}.
#' @return A ggplot object displaying both target and envelope densities.
#' @export
plot.envelope_class <- function(object, width = 1, log_scale = T, lower = 0, upper = 0.5) {

  top_point <- optimise(object$log_target_density, interval = c(lower, upper), maximum = TRUE)$maximum

  gridpoints <- seq(0, top_point+width, length.out = 512)

  transform <- function(x){
    if (!log_scale) {
      return(exp(x))
    }
    x
  }

  df <- data.frame(
    gridpoints = gridpoints,
    target = object$log_target_density(gridpoints),
    envelope = object$log_envelope_density(gridpoints)
  ) %>% pivot_longer(cols = -gridpoints, names_to = "Description", values_to = "y")

  # Plot the function
  ggplot(df, aes(x = gridpoints, y = transform(y), color = Description)) +
    geom_line() +
    labs(title = "Plot of density",subtitle = paste0("Estimated acceptance rate is ",round(object$est_alpha*100, 2),"%"), x = "y", y = "") +
    scale_y_continuous(labels = ifelse(!log_scale, label_scientific(digits = 2),label_number(accuracy = 1)))
}


plot_diff <- function(object, ...) {
  UseMethod("plot_diff")  # Tells R to dispatch methods based on the object's class
}


#' Plot Difference Between Target and Envelope Densities
#'
#' Plots the difference between the target and envelope densities for an envelope class object.
#'
#' @param object An object of class \code{envelope_class} containing the target and envelope density functions.
#' @param width A numeric value defining additional range width for the plot beyond the peak point. Default is \code{1}.
#' @param log_scale A logical value for plotting on a logarithmic scale. Default is \code{TRUE}.
#' @param lower, upper Numeric values for the optimization and density evaluation bounds. Default is \code{0} and \code{0.5}.
#' @return A ggplot object showing the difference between target and envelope densities.
#' @export
plot_diff.envelope_class <- function(object, width = 1, log_scale = T, lower = 0, upper = 0.5) {

  top_point <- optimise(object$log_target_density, interval = c(lower, upper), maximum = TRUE)$maximum

  gridpoints <- seq(0, top_point+width, length.out = 512)

  transform <- function(x){
    if (!log_scale) {
      return(exp(x))
    }
    x
  }

  df <- data.frame(
    gridpoints = gridpoints,
    target = object$log_target_density(gridpoints),
    envelope = object$log_envelope_density(gridpoints)
  ) %>%
    mutate(y = transform(envelope)-transform(target))

  # Plot the function
  ggplot(df, aes(x = gridpoints, y = y)) +
    geom_line() +
    scale_y_continuous(labels = ifelse(!log_scale, label_scientific(digits = 2),label_number(accuracy = 1))) +
    labs(title = "Difference between target and envelope", x = "x", y = "envelope-target")
}


#' Plot Rejection Sampling Results
#'
#' Visualizes samples from a target density, with an optional normalized density curve.
#'
#' @param object An object of class \code{rejection_sampler} containing sampling results and density functions.
#' @param density_function An optional function for adding a density curve. Default is \code{NULL}.
#' @param add_normalized_density A logical value to include a normalized target density curve. Default is \code{TRUE}.
#' @param lower, upper Numeric values defining the range for integrating the target density. Default is \code{0} and \code{Inf}.
#' @return A ggplot object showing sampled density with optional density curves.
#' @export
plot.rejection_sampler <- function(object, density_function = NULL, add_normalized_density = T, lower = 0, upper = Inf) {

  plot <- ggplot() +
    labs(x = "y", y = "Density") +
    theme_bw()

  # Add histogram if requested
  plot <- plot + geom_histogram(aes(x = object$samples, y = after_stat(density)),
                                color = "white", fill = "steelblue", bins = 30)

  # Capture the modified plot with density line
  if(add_normalized_density){

    target_density <- function(x) exp(object$log_target_density(x))

    int_val <- integrate(target_density, lower, upper)

    normalized_density <- function(x) target_density(x)/int_val$value

    x <- seq( min(object$samples)*ifelse(min(object$samples)>0, 0.9, 1.1), max(object$samples)*ifelse(max(object$samples)>0, 0.9, 1.1), length.out = 512)

    plot <- plot +
      geom_line(aes(x = x, y = normalized_density(x)), linewidth = 1.5, color = "red")
  }


  if(!is.null(density_function)) {
    x <- seq( min(object$samples)*ifelse(min(object$samples)>0, 0.9, 1.1), max(object$samples)*ifelse(max(object$samples)>0, 0.9, 1.1), length.out = 512)


    plot <- plot +
      geom_line(aes(x = x, y = density_function(x)), linewidth = 1.5, color = "green")
  }





  return(plot +
           theme(plot.background = element_rect(fill = "#fafafa"),
                 panel.border = element_blank()) +
           labs(title = "Samples from target with normalized density",
                subtitle = paste0("Acceptance rate is ", round(object$acceptance_rate*100, 2),"%") ))
}


#' Plot Combined Envelope and Sampling Results
#'
#' Creates a combined plot showing the envelope and rejection sampling results, using a wrapper class.
#'
#' @param object An object of class \code{wrapper_class} containing envelope and sampling objects.
#' @param log_scale A logical value to determine if the plot should use a logarithmic scale. Default is \code{FALSE}.
#' @param add_normalized_density A logical value to add a normalized density curve. Default is \code{TRUE}.
#' @param ... Additional arguments passed to the plotting functions.
#' @return A ggplot object showing a combined plot of envelope and sampling results.
#' @export
plot.wrapper_class <- function(object, log_scale = F, add_normalized_density = T, ...) {
  plot_envelope <- plot(wraps$envelope_class, log_scale = log_scale, ...)
  plot_samples <- plot(wraps$rejection_sampler, add_normalized_density)

  plot_envelope /
    plot_samples
}
