% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PBM.R
\name{pbm.gextract_local}
\alias{pbm.gextract_local}
\title{Extract local PBM scores / response}
\usage{
pbm.gextract_local(pbm, intervals, response = FALSE, scaling_q = 0.995)
}
\arguments{
\item{pbm}{A PBM object}

\item{intervals}{Genomic intervals to extract}

\item{response}{Logical, whether to compute the response (default is FALSE)}

\item{scaling_q}{Scaling quantile for normalization (default is NULL)}
}
\value{
A matrix of extracted local energies / response
}
\description{
This function extracts local PBM scores for genomic intervals and optionally computes the response.
}
