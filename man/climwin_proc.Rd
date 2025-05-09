% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/climwin_proc.R
\name{climwin_proc}
\alias{climwin_proc}
\title{Conduct climate window analysis}
\usage{
climwin_proc(
  biol_data,
  clim_data,
  ID,
  randwin = FALSE,
  seednum = 1302,
  repeats = 200,
  plot_check = FALSE,
  cinterval = "week",
  out_clim = NULL,
  stat = "mean",
  oneGrid = TRUE,
  explanYear = TRUE,
  startWindow = 0,
  endWindow = 52,
  RefMon = NA,
  weatherVar = NA
)
}
\arguments{
\item{biol_data}{Data frame with trait data for a given population and species.}

\item{clim_data}{Either raster stack of climatic data for the region that encompasses
the study location or a dataframe with the local weather from the closest weather station.}

\item{ID}{Numeric giving a unique ID of the current dataset for a given population and species.}

\item{randwin}{Logical (TRUE/FALSE). Should \code{\link[climwin]{randwin}}
be run together with \code{\link[climwin]{slidingwin}}?}

\item{seednum}{Numeric specifying the seed number, for reproducibility.}

\item{repeats}{If randwin is run, the number of times that data will be randomised.}

\item{plot_check}{Logical (TRUE/FALSE). Whether to display
one year (the first in the series) of the climatic data
 together with the study location.}

\item{cinterval}{A character specifying the resolution at which climate window
analysis will be conducted, defaults to 'week' (weekly resolution). For more details
see the same option in the function \code{\link[climwin]{slidingwin}}.}

\item{out_clim}{Character specifying the library on the path where
to save the results of climate window analysis}

\item{stat}{Character specifying which statistics to use for
aggregating the climatic data (see also the same option in
\code{\link[climwin]{slidingwin}}).}

\item{oneGrid}{Logical (TRUE/FALSE). Whether to extact climatic data from
a single grid cell where into which the study location falls or to use mean of
five cells: the focal one and four neighbours.}

\item{explanYear}{Logical (TRUE/FALSE). Whether to include year as an
explanatory variable in the baseline formula for window analysis.}

\item{startWindow}{A numeric specifying the closest number of time intervals
(set by cinterval) back from the reference day. See also the option range in
the function \code{\link[climwin]{slidingwin}}.}

\item{endWindow}{A numeric specifying the furthest number of time intervals
(set by cinterval) back from the reference day. See also the option range in
the function \code{\link[climwin]{slidingwin}}.}

\item{RefMon}{A numeric specifying the month for the absolute window
to be used for all species. If NA is specified then a species-specific window is used for
each dataset using the available data on phenology and the period when morphology
was recorded.}

\item{weatherVar}{Character corresponding to the name of the climatic variable
in the data frame provided as clim_data.}
}
\value{
A tibble with six columns if randwin is set to FALSE
and seven columns if randwin is set TRUE. The columns are: "ID"
- study id, "Species" - study species, "Trait" - the trait considered,
"climwin_output" - a list returned by slidingwin(),
"randwin_output" - a list returned by randwin (if set to TRUE),
"clim_data" - a data frame with climate values,
"biol_data" - a data frame with traits, demographic rates
and population size.
}
\description{
\code{climwin_proc} performs climate window analysis for
specified traits and climate data and returns the best
identified window together with results of randomisation
tests (if asked for)
}
\examples{
# ATTENTION: DO NOT RUN! takes long time
\dontrun{
biol_noSea <- prep_subset(data = data_biol, Seabird = FALSE)
# keep only EU countries
biol_eu <- droplevels(subset(biol_noSea$subdata[[1]],
! Country \%in\% c('Antarctica', 'Australia',
                'Canada', 'Falkland Islands',
                'Greenland', 'Mexico',
                'New Zealand', 'South Africa',
                'South Atlantic Ocean',
                'South Georgia', 'Svalbard',
                'Taiwan', 'USA', 'Venezuela')))
meanT <- raster::stack(x = system.file("extdata",
"tg_ens_mean_0.1deg_reg_v18.0e.nc", package="sTraitChange"))
test_rand <- climwin_proc(biol_data = biol_eu,
                          clim_data = meanT, ID = 1,
                          randwin = FALSE, seednum = 1302,
                          repeats = 30, plot_check = FALSE,
                          out_clim = tempdir(), # attention: for this example we
                                     # write the data to a temporary directory,
                                     # to check its location type tempdir()
                          cinterval = 'month',
                          stat = 'mean',
                          startWindow = 0, endWindow = 12,
                          oneGrid = FALSE, explanYear = TRUE,
                          RefMon = NA, weatherVar = NA)}
message('Temporary directory is located at', tempdir())
}
