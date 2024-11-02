#' Extract results of the fitted structural equation model
#'
#' \code{get_res_SEM} extracts results of the fitted
#'  structural equation model
#'
#' @param mod_obj fitted psem object, returned with function
#' \code{\link{fit_mod}}
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @return A tibble with 29 columns that represent the metadata for a
#' given study ID and the results of SEM: Cstat, dTable, R2 of all
#' equations and path coefficients.
#'
#' @examples
#'
#'  # prepare the data to fit the model to one study only
#'  subs <- droplevels(dataSEM[dataSEM$ID == 7, ])
#'  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1))
#'  consec_yrs <- merge(full_NA, subs, by = 'Year', all= TRUE)
#'
#'  data_GR <- consec_yrs %>%
#'  dplyr::mutate(Pop_mean_lag = c(.data$Pop_mean[-1], NA)) %>%
#'  dplyr::mutate(GR = log(.data$Pop_mean_lag / .data$Pop_mean)) %>%
#'  dplyr::filter(!is.na(.data$GR) & !is.na(.data$Trait_mean) &
#'                !is.na(.data$Demog_rate_mean) & !is.na(.data$Pop_mean)) %>%
#'  dplyr::mutate(det_Clim = stats::resid(stats::lm(Clim ~ Year,
#'                data = .))) %>%
#'  dplyr::mutate(across(where(is.array), as.numeric))
#'
#' # fit the model
#' test <- fit_mod(biol_data = data_GR, ID = 7,
#'                 DD = 'n_effectGR', weight = TRUE,
#'                 correlation = TRUE,
#'                 Trait = TRUE,
#'                 simpleSEM = TRUE)
#'  # extract results
#' res_SEM <- get_res_SEM(test)

get_res_SEM <- function(mod_obj){

  suppressWarnings(summary_mod <- summary(mod_obj, .progressBar = F))
  dTable <- summary_mod$dTable
  Cstat <- summary_mod$Cstat
  R2_Relation <- piecewiseSEM::rsquared(mod_obj)
  coefs <- piecewiseSEM::coefs(mod_obj)

  AIC <- stats::AIC(mod_obj, aicc = TRUE)
  data <- mod_obj$data

  # add timeseries duration
  data <- data %>%
    dplyr::group_by(data$ID) %>%
    dplyr::mutate(NYears = dplyr::n()) %>%
    dplyr::ungroup()

  res <- tibble::tibble(ID = unique(data$ID),
                        Study_Authors = unique(data$Study_Authors),
                        Species = unique(data$Species),
                        Location = unique(data$Location),
                        Country = unique(data$Country),
                        Continent = unique(data$Continent),
                        Longitude = unique(data$Longitude),
                        Latitude = unique(data$Latitude),
                        Taxon = unique(data$Taxon),
                        BirdType = unique(data$BirdType),
                        Trait_Categ = unique(data$Trait_Categ),
                        Trait = unique(data$Trait),
                        Demog_rate_Categ = unique(data$Demog_rate_Categ),
                        Demog_rate = unique(data$Demog_rate),
                        Count = unique(data$Count),
                        Nyears = unique(data$NYears),
                        WinDur = unique(data$WinDur),
                        Ref.day = unique(data$Ref.day),
                        Ref.month = unique(data$Ref.month),
                        WindowClose = unique(data$WindowClose),
                        deltaAIC = unique(data$deltaAIC),
                        Pvalue = unique(data$Pvalue),
                        Trait_ageClass = unique(data$Trait_ageClass),
                        dTable = list(dTable),
                        Cstat = list(Cstat),
                        R2 = list(R2_Relation),
                        coefs = list(coefs),
                        AIC = AIC,
                        WeathQ = unique(data$Weather_data_quality))

  return(res)
}
