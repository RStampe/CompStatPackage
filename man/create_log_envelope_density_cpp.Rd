% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{create_log_envelope_density_cpp}
\alias{create_log_envelope_density_cpp}
\title{Log Envelope Density Calculation (C++)}
\usage{
create_log_envelope_density_cpp(x, z, a, b)
}
\arguments{
\item{x}{A numeric vector of values for which to calculate the envelope density.}

\item{z}{A numeric vector for interval boundaries.}

\item{a}{A numeric vector of coefficients.}

\item{b}{A numeric vector of coefficients.}
}
\value{
A numeric vector of log envelope density values.
}
\description{
Calculates a piecewise linear log envelope density using input parameters.
}
