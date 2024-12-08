% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM.R
\name{calculate_mle}
\alias{calculate_mle}
\title{Maximum Likelihood Estimation}
\usage{
calculate_mle(x, w, nu)
}
\arguments{
\item{x}{A numeric vector of data.}

\item{w}{A numeric vector of weights.}

\item{nu}{Degrees of freedom parameter.}
}
\value{
A numeric vector with estimated \code{mu} and \code{sigma2}.
}
\description{
Calculates the maximum likelihood estimators for \code{mu} and \code{sigma2}.
}
