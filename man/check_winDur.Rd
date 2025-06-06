% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_winDur.R
\name{check_winDur}
\alias{check_winDur}
\title{Check whether the window duration is within realistic range}
\usage{
check_winDur(climwin_out, clim, biol_data, MinDur, MaxDur, deltaThresh = -7)
}
\arguments{
\item{climwin_out}{A list containing the results of climatic window analyses,
returned by \code{\link[climwin]{slidingwin}} called within \code{\link{climwin_proc}}.}

\item{clim}{A dataframe consisting of daily climatic values,
corresponding to the clim_data element as returned by \code{\link{climwin_proc}}.}

\item{biol_data}{Data frame with trait data for a given
population and species.}

\item{MinDur}{A numeric specifying the minimum allowed window duration.}

\item{MaxDur}{A numeric specifying the maximum allowed window duration.}

\item{deltaThresh}{A numeric specifying the minimum deltaAIC value for which the
model will still be judged acceptable for selecting the window duration.}
}
\value{
If the climatic window within the specified range is found and the
respective model has deltaAIC at least -7 then the function returns a
dataframe consisting of biological data and an additional column 'Clim'
that for each year specifies the climatic values obtained according to
the selected model. If either the climatic window satisfying the specified range
is not found or / and the respective model has deltaAIC > -7, the function returns
NULL.
}
\description{
\code{check_winDur} checks whether the window duration of the best
selected model is within the specified range and if not, finds the
model with such a window (given that the model with such window
is still better than the one without climatic variable as explanatory)
}
