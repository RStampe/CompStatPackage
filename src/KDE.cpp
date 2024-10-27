#include <Rcpp.h>
#include <unistd.h>
#include <cmath>
using namespace Rcpp;

double epanechnikov_kernel_single_cpp(double x) {
  if (-1 < x && x < 1)
    return (1 - x*x) * 0.75;
  else
    return 0;
}

//' Epanechnikov Kernel Function (C++)
//'
//' Calculates the Epanechnikov kernel for a given vector of input values.
//'
//' @param x A numeric vector of values for which to calculate the Epanechnikov kernel.
//' @return A numeric vector of kernel values, with zeros for inputs outside the range -1, 1.
//' @export
// [[Rcpp::export]]
NumericVector epanechnikov_kernel_cpp(const NumericVector& x) {
  // if problem with int use unsigned long long
  int n = x.size();

  NumericVector result(n);
  for(int i = 0; i < n; ++i) {
    result[i] = epanechnikov_kernel_single_cpp(x[i]);
  }
  return result;
}

double gaussian_kernel_single_cpp(double x) {
  const double inv_sqrt_2pi = 0.3989422804014337; // Precompute 1 / sqrt(2 * PI)
  return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

//' Gaussian Kernel Function (C++)
//'
//' Calculates the Gaussian (Normal) kernel for a given vector of input values.
//'
//' @param x A numeric vector of values for which to calculate the Gaussian kernel.
//' @return A numeric vector of Gaussian kernel values for each element in `x`.
//' @export
// [[Rcpp::export]]
NumericVector gaussian_kernel_cpp(NumericVector x) {
  // if problem with int use unsigned long long
  int n = x.size();

  NumericVector result(n);
  for(int i = 0; i < n; ++i) {
    result[i] = gaussian_kernel_single_cpp(x[i]);
  }
  return result;
}

//' Epanechnikov Kernel Density Estimation (C++)
//'
//' Computes kernel density estimates at specified grid points using the Epanechnikov kernel.
//' The function is parallelized with OpenMP for improved performance.
//'
//' @param grid_points A numeric vector of grid points where density estimates are evaluated.
//' @param vec_obs A numeric vector of observed data points.
//' @param bandwidth A positive numeric value specifying the smoothing bandwidth.
//' @return A list with two components:
//'   - `x`: Grid points where the density is estimated.
//'   - `y`: Corresponding density estimates using the Epanechnikov kernel.
//' @export
// [[Rcpp::export]]
List kernel_cpp_epanechnikov(const NumericVector& grid_points, const NumericVector& vec_obs, double bandwidth) {
  int n_grid = grid_points.size();
  int n_obs = vec_obs.size();

  NumericVector vec_density_values(n_grid);

  // Precompute 1.0 / (bandwidth * n_obs) to avoid redundant division
  double inv_bandwidth_n_obs = 1.0 / (bandwidth * n_obs);

  for (int i = 0; i < n_grid; i++) {
    double density_value = 0.0;

    // Loop over each observation
    for (int j = 0; j < n_obs; j++) {
      double diff = (grid_points[i] - vec_obs[j]) / bandwidth;
      density_value += epanechnikov_kernel_single_cpp(diff);
    }

    // Compute the final density value
    vec_density_values[i] = density_value * inv_bandwidth_n_obs;
  }

  // Return the result as a list with grid points and density values
  return List::create(Named("x") = grid_points,
                      Named("y") = vec_density_values);
}


//' Gaussian Kernel Density Estimation (C++)
//'
//' Computes kernel density estimates at specified grid points using the Gaussian kernel.
//' The function is parallelized with OpenMP for optimized performance on large data sets.
//'
//' @param grid_points A numeric vector of grid points where density estimates are evaluated.
//' @param vec_obs A numeric vector of observed data points.
//' @param bandwidth A positive numeric value specifying the smoothing bandwidth.
//' @return A list with two components:
//'   - `x`: Grid points where the density is estimated.
//'   - `y`: Corresponding density estimates using the Gaussian kernel.
//' @export
// [[Rcpp::export]]
List kernel_cpp_gaussian(const NumericVector& grid_points, const NumericVector& vec_obs, double bandwidth) {
  int n_grid = grid_points.size();
  int n_obs = vec_obs.size();

  NumericVector vec_density_values(n_grid);

  // Precompute 1.0 / (bandwidth * n_obs) to avoid redundant division
  double inv_bandwidth_n_obs = 1.0 / (bandwidth * n_obs);

  // Loop over each grid point with OpenMP
  for (int i = 0; i < n_grid; i++) {
    double density_value = 0.0;

    // Loop over each observation
    for (int j = 0; j < n_obs; j++) {
      double diff = (grid_points[i] - vec_obs[j]) / bandwidth;
      density_value += gaussian_kernel_single_cpp(diff);
    }

    // Compute the final density value
    vec_density_values[i] = density_value * inv_bandwidth_n_obs;
  }

  // Return the result as a list with grid points and density values
  return List::create(Named("x") = grid_points,
                      Named("y") = vec_density_values);
}


//' Data Binning Function (C++)
//'
//' Bins observed data points within a specified range, dividing the range into equal-width bins.
//'
//' @param vec_obs A numeric vector of observed data points.
//' @param lower The lower bound of the binning range.
//' @param upper The upper bound of the binning range.
//' @param number_of_bins An integer specifying the number of bins.
//' @return A numeric vector with counts of observations in each bin.
//' @export
// [[Rcpp::export]]
NumericVector binning_cpp(const NumericVector& vec_obs, double lower, double upper, int number_of_bins) {
  double delta = (upper - lower) / (number_of_bins - 1);
  std::vector<double> binned_obs(number_of_bins, 0.0);

  // Parallel loop over vec_obs to accumulate counts
  for (int j = 0; j < vec_obs.size(); j++) {
    int i = std::floor((vec_obs[j] - lower) / delta + 0.5);

    // Ensure bin index is within bounds
    if (i >= 0 && i < number_of_bins) {
      binned_obs[i] += 1;
    }
  }

  // Sum all bin counts
  double total_count = 0.0;
  for (int i = 0; i < number_of_bins; i++) {
    total_count += binned_obs[i];
  }

  // Normalize bin counts
  NumericVector normalized_bins(number_of_bins);
  if (total_count > 0) {
    for (int i = 0; i < number_of_bins; i++) {
      normalized_bins[i] = binned_obs[i] / total_count;
    }
  }

  return normalized_bins;
}

//' Epanechnikov Kernel Density Estimation with Binning (C++)
//'
//' Computes kernel density estimates at specified grid points using the Epanechnikov kernel
//' with a binning approach for improved performance.
//'
//' @param grid_points A numeric vector of grid points where density estimates are evaluated.
//' @param vec_obs A numeric vector of observed data points.
//' @param bandwidth A positive numeric value specifying the smoothing bandwidth.
//' @return A list with two components:
//'   - `x`: Grid points where the density is estimated.
//'   - `y`: Corresponding density estimates using the Epanechnikov kernel and binning.
//' @export
// [[Rcpp::export]]
List kernel_cpp_epanechnikov_binning(const NumericVector& grid_points,
                                     const NumericVector& vec_obs,
                                     double bandwidth) {
  int i, j;
  int grid_length = grid_points.length();
  double grid_lower = grid_points[0];
  double grid_upper = grid_points[grid_length - 1];
  double grid_diff = grid_points[1] - grid_points[0];
  NumericVector y(grid_length);

  // Corrected: Pass grid_upper to binning_cpp
  NumericVector weights = binning_cpp(vec_obs, grid_lower, grid_upper, grid_length);
  int max_nonzero_i = floor(bandwidth / grid_diff);

  // Pad weights with zeros at the beginning and end
  NumericVector weights_padded(grid_length + 2 * max_nonzero_i);
  for (i = 0; i < grid_length; i++) {
    weights_padded[max_nonzero_i + i] = weights[i];
  }

  // Compute kernel_evals
  NumericVector kernel_evals(max_nonzero_i + 1);
  kernel_evals[0] = 0.75; // Epanechnikov kernel at 0
  for (i = 1; i <= max_nonzero_i; i++) {
    kernel_evals[i] = epanechnikov_kernel_single_cpp((i * grid_diff) / bandwidth);
  }

  // Construct kernel_vec by concatenating reversed kernel_evals[-1] and kernel_evals
  int total_kernel_length = 2 * max_nonzero_i + 1;
  NumericVector kernel_vec(total_kernel_length);
  for (i = 0; i < max_nonzero_i; i++) {
    kernel_vec[i] = kernel_evals[max_nonzero_i - i];
  }
  kernel_vec[max_nonzero_i] = kernel_evals[0];
  for (i = 1; i <= max_nonzero_i; i++) {
    kernel_vec[max_nonzero_i + i] = kernel_evals[i];
  }

  // Perform convolution
  for (i = 0; i < grid_length; i++) {
    double sum = 0.0;
    for (j = 0; j < total_kernel_length; j++) {
      sum += weights_padded[i + j] * kernel_vec[j];
    }
    y[i] = sum / bandwidth;
  }

  // Return the result as a list with grid points and density values
  return List::create(Named("x") = grid_points,
                      Named("y") = y);
}



//' Gaussian Kernel Density Estimation with Binning (C++)
//'
//' Computes kernel density estimates at specified grid points using the Gaussian kernel
//' with a binning approach for improved performance.
//'
//' @param grid_points A numeric vector of grid points where density estimates are evaluated.
//' @param vec_obs A numeric vector of observed data points.
//' @param bandwidth A positive numeric value specifying the smoothing bandwidth.
//' @return A list with two components:
//'   - `x`: Grid points where the density is estimated.
//'   - `y`: Corresponding density estimates using the Gaussian kernel and binning.
//' @export
// [[Rcpp::export]]
List kernel_cpp_gaussian_binning(const NumericVector& grid_points,
                                 const NumericVector& vec_obs,
                                 double bandwidth) {
  int i, j;
  int grid_length = grid_points.length();
  double grid_lower = grid_points[0];
  double grid_upper = grid_points[grid_length - 1];
  double grid_diff = grid_points[1] - grid_points[0];
  NumericVector y(grid_length);

  // Corrected: Pass grid_upper to binning_cpp
  NumericVector weights = binning_cpp(vec_obs, grid_lower, grid_upper, grid_length);
  int max_nonzero_i = floor(bandwidth / grid_diff);

  // Pad weights with zeros at the beginning and end
  NumericVector weights_padded(grid_length + 2 * max_nonzero_i);
  for (i = 0; i < grid_length; i++) {
    weights_padded[max_nonzero_i + i] = weights[i];
  }

  // Compute kernel_evals
  NumericVector kernel_evals(max_nonzero_i + 1);
  kernel_evals[0] = 0.75; // Epanechnikov kernel at 0
  for (i = 1; i <= max_nonzero_i; i++) {
    kernel_evals[i] = gaussian_kernel_single_cpp((i * grid_diff) / bandwidth);
  }

  // Construct kernel_vec by concatenating reversed kernel_evals[-1] and kernel_evals
  int total_kernel_length = 2 * max_nonzero_i + 1;
  NumericVector kernel_vec(total_kernel_length);
  for (i = 0; i < max_nonzero_i; i++) {
    kernel_vec[i] = kernel_evals[max_nonzero_i - i];
  }
  kernel_vec[max_nonzero_i] = kernel_evals[0];
  for (i = 1; i <= max_nonzero_i; i++) {
    kernel_vec[max_nonzero_i + i] = kernel_evals[i];
  }

  // Perform convolution
  for (i = 0; i < grid_length; i++) {
    double sum = 0.0;
    for (j = 0; j < total_kernel_length; j++) {
      sum += weights_padded[i + j] * kernel_vec[j];
    }
    y[i] = sum / bandwidth;
  }

  // Return the result as a list with grid points and density values
  return List::create(Named("x") = grid_points,
                      Named("y") = y);
}





