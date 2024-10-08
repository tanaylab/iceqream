% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IQFeatureGroup.R
\docType{class}
\name{IQFeatureGroup-class}
\alias{IQFeatureGroup-class}
\alias{IQFeatureGroup}
\alias{initialize,IQFeatureGroup-method}
\alias{show,IQFeatureGroup-method}
\title{IQFeatureGroup class}
\usage{
\S4method{initialize}{IQFeatureGroup}(.Object, features, compute_func, size)

\S4method{show}{IQFeatureGroup}(object)
}
\arguments{
\item{.Object}{The IQFeatureGroup object to initialize.}

\item{features}{A list of IQSeqFeature objects.}

\item{compute_func}{A function to compute all features in the group at once.}

\item{size}{The size of the sequences to use for feature computation.}

\item{object}{An IQFeatureGroup object}
}
\value{
An initialized IQFeatureGroup object.
}
\description{
This class represents a group of related IQ features, such as all dinucleotide features.

This method defines how an IQFeatureGroup object should be displayed.
}
\section{Slots}{

\describe{
\item{\code{features}}{A list of IQSeqFeature objects.}

\item{\code{compute_func}}{A function to compute all features in the group at once.}

\item{\code{size}}{The size of the sequences to use for feature computation (numeric).}
}}

