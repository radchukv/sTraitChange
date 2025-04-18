% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_lab_name.R
\name{gen_specific_labs}
\alias{gen_specific_labs}
\title{Turn general axes labels into specific ones}
\usage{
gen_specific_labs(
  vect,
  Clim = "Temperature",
  Demog_rate = "Survival",
  Trait_categ = "Phenology"
)
}
\arguments{
\item{vect}{character specifying the general axis label that is to be converted
into a specific one.}

\item{Clim}{Character specifying the level of the climatic variable, for which analyses
were conducted.}

\item{Demog_rate}{Character specifying the level of the demographic rate for which
analyses were conducted.}

\item{Trait_categ}{Character specifying the level of the trait, for which analyses were
conducted.}
}
\value{
A vector of length one with the axis label.
}
\description{
\code{gen_specific_labs} prepares the axes labels for the forest
and 'relation' plots. It replaces general names in axes (e.g. demographic rate)
by specific ones (e.g. survival). This is an internal function
used by \code{\link{plot_forest}}.
}
