% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_res_SEM.R
\name{get_res_SEM}
\alias{get_res_SEM}
\title{Extract results of the fitted structural equation model}
\usage{
get_res_SEM(mod_obj)
}
\arguments{
\item{mod_obj}{fitted psem object, returned with function
\code{\link{fit_mod}}}
}
\value{
A tibble with 29 columns that represent the metadata for a
given study ID and the results of SEM: Cstat, dTable, R2 of all
equations and path coefficients.
}
\description{
\code{get_res_SEM} extracts results of the fitted
 structural equation model
}
\examples{

 # prepare the data to fit the model to one study only
 subs <- droplevels(dataSEM[dataSEM$ID == 7, ])
 full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1))
 consec_yrs <- merge(full_NA, subs, by = 'Year', all= TRUE)

 data_GR <- consec_yrs \%>\%
 dplyr::mutate(Pop_mean_lag = c(.data$Pop_mean[-1], NA)) \%>\%
 dplyr::mutate(GR = log(.data$Pop_mean_lag / .data$Pop_mean)) \%>\%
 dplyr::filter(!is.na(.data$GR) & !is.na(.data$Trait_mean) &
               !is.na(.data$Demog_rate_mean) & !is.na(.data$Pop_mean)) \%>\%
 dplyr::mutate(det_Clim = stats::resid(stats::lm(Clim ~ Year,
               data = .))) \%>\%
 dplyr::mutate(across(where(is.array), as.numeric))

# fit the model
test <- fit_mod(biol_data = data_GR, ID = 7,
                DD = 'n_effectGR', weight = TRUE,
                correlation = TRUE,
                Trait = TRUE,
                simpleSEM = TRUE)
 # extract results
res_SEM <- get_res_SEM(test$models_list[[1]])
}
