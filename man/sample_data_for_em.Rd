% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EM.R
\name{sample_data_for_em}
\alias{sample_data_for_em}
\title{Generate a sample of data with parameters}
\usage{
sample_data_for_em(n, parameters = c(0.5, 1, 1))
}
\arguments{
\item{n}{An integer specifying the number of samples to generate.}

\item{parameters}{A numeric vector of length three for \code{mu}, \code{sigma2}, and \code{nu}.}
}
\value{
A list with elements:
\item{x}{A numeric vector of sampled data.}
\item{w}{A numeric vector of chi-squared distributed values.}
}
\description{
This function generates a random sample of data \code{x} from a normal distribution
with mean \code{mu}, variance \code{sigma2}, and degrees of freedom \code{nu}. It also provides
a \code{w} vector of chi-squared distributed values.
}
