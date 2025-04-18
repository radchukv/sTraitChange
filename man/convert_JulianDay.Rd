% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_JulianDay.R
\name{convert_JulianDay}
\alias{convert_JulianDay}
\title{Convert phenological traits that are not specified
in Julian days to Julian days}
\usage{
convert_JulianDay(biol_data)
}
\arguments{
\item{biol_data}{Data frame with trait data for a given population and species.}
}
\value{
A data frame with columns as those in the input data
(biol_data), but with any phenological traits that were specified
not in Julian Days converted to Julian Days.
}
\description{
\code{convert_JulianDay} performs standardization of phenological
traits by converting those that are specified from other
starting day than the 1st of January to Julian days
}
\examples{
out_biol <- convert_JulianDay(biol_data = data_biol)

}
