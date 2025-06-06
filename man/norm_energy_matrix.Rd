% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/energy-utils.R
\name{norm_energy_matrix}
\alias{norm_energy_matrix}
\title{Normalize Energy Matrix}
\usage{
norm_energy_matrix(
  x,
  dataset_x = x,
  min_energy = -7,
  q = 1,
  norm_energy_max = 10
)
}
\arguments{
\item{x}{The input matrix to be normalized.}

\item{dataset_x}{The reference dataset matrix used for normalization.}

\item{min_energy}{The minimum energy value to be assigned after normalization. Default is -7.}

\item{q}{The quantile value used for calculating the maximum value in the reference dataset. Default is 1.}

\item{norm_energy_max}{The maximum value to which the normalized energy values are scaled. Default is 10.}
}
\value{
A normalized energy matrix with the same dimensions as the input matrix.
}
\description{
This function normalizes an energy matrix by applying logarithmic transformation and scaling.
}
\examples{
# Example usage:
data <- matrix(rnorm(100), nrow = 10)
normalized_data <- norm_energy_matrix(data, data, min_energy = -7, q = 1, norm_energy_max = 10)

}
