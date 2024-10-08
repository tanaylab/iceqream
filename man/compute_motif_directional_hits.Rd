% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial.R
\name{compute_motif_directional_hits}
\alias{compute_motif_directional_hits}
\title{Compute motif directional hits}
\usage{
compute_motif_directional_hits(pssm, intervals, size, pwm_threshold = 7)
}
\arguments{
\item{pssm}{A position-specific scoring matrix (PSSM) object.}

\item{intervals}{A data frame containing the genomic intervals to search for motif hits.}

\item{size}{The size of the intervals.}

\item{pwm_threshold}{The threshold for the PSSM score to consider a hit.}
}
\value{
A data frame containing the number of hits of the motif in each interval, separated by strand and direction. Column names are 'ml' (minus, left), 'mr' (minus, right), 'pl' (plus, left), 'pr' (plus, right).
}
\description{
This function computes the number of hits of a given motif in a set of intervals, separated by strand and direction.
}
