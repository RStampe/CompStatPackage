#include <Rcpp.h>
#include <cmath>  // For exp()

using namespace Rcpp;

// ==========================================
// ========== Initializing Target ===========
// ==========================================

// Declare global vectors
NumericVector x, z;

//' @name set_data_vectors
 //' @title Set Data Vectors (C++)
 //'
 //' Initializes the global vectors `x` and `z` with input values from R.
 //'
 //' @param x_in A numeric vector to set the global vector `x`.
 //' @param z_in A numeric vector to set the global vector `z`.
 //' @return None. This function is used for its side effects.
 //' @export
 // [[Rcpp::export]]
 void set_data_vectors(NumericVector x_in, NumericVector z_in) {
   x = x_in;
   z = z_in;
 }

//' @name log_target_density_single_obs
 //' @title Log Target Density for Single Observation (C++)
 //'
 //' Computes the log target density for a single observation, given a numeric value.
 //'
 //' @param y A numeric value representing a single observation.
 //' @return The log target density for \code{y}. Returns \code{-Inf} if \code{y < 0}.
 //' @export
 double log_target_density_single_obs(double y) {
   if (y < 0) {
     return -std::numeric_limits<double>::infinity();
   }
   int m = 100;
   double result = 0.0;
   for (int j = 0; j < m; j++) {
     double term = y * z[j] * x[j] - std::exp(y * x[j]);
     result += term;
   }
   return result;
 }

//' @name log_target_density_cpp
 //' @title Log Target Density for Vector of Observations (C++)
 //'
 //' Calculates the log target density for each value in a numeric vector.
 //'
 //' @param y A numeric vector of values for which to calculate the log target density.
 //' @return A numeric vector of log target density values.
 //' @export
 // [[Rcpp::export]]
 NumericVector log_target_density_cpp(NumericVector y) {
   int n = y.size();
   NumericVector result(n);
   for (int i = 0; i < n; i++) {
     result[i] = log_target_density_single_obs(y[i]);
   }
   return result;
 }

//' @name d_log_target_density_single_obs
 //' @title Derivative of Log Target Density for Single Observation (C++)
 //'
 //' Computes the derivative of the log target density for a single observation.
 //'
 //' @param y A numeric value representing a single observation.
 //' @return The derivative of the log target density for \code{y}. Returns \code{-Inf} if \code{y < 0}.
 //' @export
 double d_log_target_density_single_obs(double y) {
   if (y < 0) {
     return -std::numeric_limits<double>::infinity();
   }
   int m = x.size();
   double result = 0.0;
   for (int j = 0; j < m; j++) {
     result += z[j] * x[j] - x[j] * std::exp(y * x[j]);
   }
   return result;
 }

//' @name d_log_target_density_cpp
 //' @title Derivative of Log Target Density for Vector of Observations (C++)
 //'
 //' Calculates the derivative of the log target density for each value in a numeric vector.
 //'
 //' @param y A numeric vector of values for which to calculate the derivative of the log target density.
 //' @return A numeric vector of derivative values of the log target density.
 //' @export
 // [[Rcpp::export]]
 NumericVector d_log_target_density_cpp(NumericVector y) {
   int n = y.size();
   NumericVector result(n);
   for (int i = 0; i < n; i++) {
     result[i] = d_log_target_density_single_obs(y[i]);
   }
   return result;
 }

// Custom findInterval function compatible with Rcpp vectors
inline int findInterval(double value, const NumericVector& breaks) {
  int idx = std::upper_bound(breaks.begin(), breaks.end(), value) - breaks.begin() - 1;
  int n = breaks.size();
  if (idx < 0) return 0;
  if (idx >= n - 1) return n - 1;
  return idx;
}


//' @name log_envelope_gaussian_density
 //' @title Log Envelope Gaussian Density (C++)
 //'
 //' Calculates the log target density for a Gaussian envelope.
 //'
 //' @param x A numeric value.
 //' @param mu Mean of the Gaussian distribution.
 //' @param sd Standard deviation of the Gaussian distribution.
 //' @param log_alpha_prime Log scaling parameter.
 //' @return The log envelope density.
 //' @export
 double log_envelope_gaussian_density(double x, double mu, double sd, double log_alpha_prime) {
   return (-((x - mu) * (x - mu)) / (2 * sd * sd) - log_alpha_prime);
 }

//' @name get_accepts_gaussian_cpp
 //' @title Gaussian Rejection Sampler Accept/Reject Decisions (C++)
 //'
 //' Determines acceptance of samples in a Gaussian rejection sampling algorithm.
 //'
 //' @param mu Mean of the Gaussian distribution.
 //' @param sd Standard deviation of the Gaussian distribution.
 //' @param log_alpha_prime Log scaling parameter for acceptance region.
 //' @param proposals A numeric vector of proposed sample values.
 //' @param u_samples A numeric vector of random uniform values for rejection sampling.
 //' @return A logical vector indicating accepted (\code{TRUE}) or rejected (\code{FALSE}) samples.
 //' @export
 // [[Rcpp::export]]
 LogicalVector get_accepts_gaussian_cpp(double mu, double sd, double log_alpha_prime, NumericVector proposals, NumericVector u_samples) {
   int n = proposals.size();
   LogicalVector accepts(n);
   for (int i = 0; i < n; ++i) {
     bool reject = u_samples[i] > std::exp(log_target_density_single_obs(proposals[i]) - log_envelope_gaussian_density(proposals[i], mu, sd, log_alpha_prime));
     accepts[i] = !reject;
   }
   return accepts;
 }

//' @name create_sampler_function_cpp
 //' @title Sampler Function with Acceptance-Rejection (C++)
 //'
 //' Samples values according to a custom acceptance-rejection method.
 //'
 //' @param n Integer specifying the number of samples.
 //' @param c Constant multiplier for scaling.
 //' @param Q A cumulative probability distribution vector.
 //' @param a A vector of coefficients.
 //' @param b A vector of coefficients.
 //' @param z A vector of control parameters.
 //' @return A numeric vector of sampled values.
 //' @export
 // [[Rcpp::export]]
 NumericVector create_sampler_function_cpp(int n, double c, const NumericVector& Q, const NumericVector& a, const NumericVector& b, const NumericVector& z) {
   NumericVector result(n);
   for (int j = 0; j < n; j++) {
     double u = R::runif(0.0, 1.0);
     double cuj = c * u;
     int i = findInterval(cuj, Q);
     result[j] = std::log(a[i] * std::exp(-b[i]) * (cuj - Q[i]) + std::exp(a[i] * z[i])) / a[i];
   }
   return result;
 }

//' @name create_log_envelope_density_cpp
 //' @title Log Envelope Density Calculation (C++)
 //'
 //' Calculates a piecewise linear log envelope density using input parameters.
 //'
 //' @param x A numeric vector of values for which to calculate the envelope density.
 //' @param z A numeric vector for interval boundaries.
 //' @param a A numeric vector of coefficients.
 //' @param b A numeric vector of coefficients.
 //' @return A numeric vector of log envelope density values.
 //' @export
 // [[Rcpp::export]]
 NumericVector create_log_envelope_density_cpp(const NumericVector& x, const NumericVector& z, const NumericVector& a, const NumericVector& b) {
   int n = x.size();
   NumericVector result(n);
   for (int j = 0; j < n; j++) {
     int i = findInterval(x[j], z);
     result[j] = a[i] * x[j] + b[i];
   }
   return result;
 }
