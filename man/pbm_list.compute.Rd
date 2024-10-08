% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PBM.R
\name{pbm_list.compute}
\alias{pbm_list.compute}
\title{Compute energies / response for a list of PBMs on given sequences}
\usage{
pbm_list.compute(
  pbm_list,
  sequences,
  response = FALSE,
  func = "logSumExp",
  normalize_energies = TRUE
)
}
\arguments{
\item{pbm_list}{A list of PBM objects}

\item{sequences}{A set of sequences on which to compute the energies.}

\item{response}{A logical flag indicating whether to compute the response. Default is FALSE.}

\item{func}{the function to use to combine the PWMs for each sequence. Either 'logSumExp' or 'max'. The default is 'logSumExp'.}

\item{normalize_energies}{A logical flag indicating whether to normalize the energies to a range of 0-10. Default is TRUE. Note that response computation requires normalized energies.}
}
\value{
A matrix containing the computed energies for each PBM
}
\description{
This function computes the energies for a list of PBMs on given sequences.
}
