% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traj-model-utils.R
\name{add_motif_models_to_traj}
\alias{add_motif_models_to_traj}
\title{Add motif models to trajectory model}
\usage{
add_motif_models_to_traj(traj_model, new_motif_models, verbose = TRUE)
}
\arguments{
\item{traj_model}{The trajectory model object to add motif models to.}

\item{new_motif_models}{A named list of motif models to add. Each element should contain a 'pssm' component and optionally a 'spat' component.}

\item{verbose}{A logical value indicating whether to display information about the R^2 after adding the motif models. Default is TRUE.}
}
\value{
The updated trajectory model object after adding the motif models.
}
\description{
This function adds specified motif models to a given trajectory model.
It updates the model, motif models, predicted difference score, model features,
coefficients, normalized energies, and features R^2 of the trajectory model.
}
