% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrajectoryModel.R
\name{predict,TrajectoryModel-method}
\alias{predict,TrajectoryModel-method}
\title{Predict TrajectoryModel on new data}
\usage{
\S4method{predict}{TrajectoryModel}(object, peak_intervals, additional_features = NULL)
}
\arguments{
\item{object}{An instance of \code{TrajectoryModel}.}

\item{peak_intervals}{data frame, indicating the genomic positions ('chrom', 'start', 'end') of each peak to predict.}

\item{additional_features}{A data frame, representing additional genomic features (e.g. CpG content, distance to TSS, etc.) for each peak. Note that NA values would be replaced with 0.}
}
\value{
A numeric vector of predicted differential accessibility scores.
}
\description{
Computes the predicted differential accessibility score between the start and end of the trajectory.
}
