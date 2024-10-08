% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interactions.R
\name{add_interactions}
\alias{add_interactions}
\title{Add interactions to a trajectory model}
\usage{
add_interactions(
  traj_model,
  interaction_threshold = 0.001,
  max_motif_n = NULL,
  max_add_n = NULL,
  lambda = 0.00001,
  alpha = 1,
  seed = 60427
)
}
\arguments{
\item{interaction_threshold}{threshold for the selecting features to create interactions. IQ learns a linear model on the features and selects the features with coefficients above this threshold. Default: 0.001}

\item{lambda}{A user supplied \code{lambda} sequence. Typical usage is to
have the program compute its own \code{lambda} sequence based on
\code{nlambda} and \code{lambda.min.ratio}. Supplying a value of
\code{lambda} overrides this. WARNING: use with care. Avoid supplying a
single value for \code{lambda} (for predictions after CV use
\code{predict()} instead).  Supply instead a decreasing sequence of
\code{lambda} values. \code{glmnet} relies on its warms starts for speed,
and its often faster to fit a whole path than compute a single fit.}

\item{alpha}{The elasticnet mixing parameter, with \eqn{0\le\alpha\le 1}.
The penalty is defined as
\deqn{(1-\alpha)/2||\beta||_2^2+\alpha||\beta||_1.} \code{alpha=1} is the
lasso penalty, and \code{alpha=0} the ridge penalty.}

\item{seed}{random seed for reproducibility.}
}
\value{
The updated trajectory model with added interactions.
}
\description{
This function adds significant interactions to a given trajectory model if they do not already exist.
It identifies significant interactions based on the provided threshold and updates the model features
with logistic features derived from these interactions. The trajectory model is then re-learned with
the new features.
}
