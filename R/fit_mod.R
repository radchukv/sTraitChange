#' Fit structural equation model with requested specifications
#'
#' \code{fit_mod} fits structural equation model with
#' requested model structure for the specified study ID
#'
#' @param biol_data Data frame with climate, trait and demographic rate
#' data for a given population and species.
#' @param ID Numeric giving a unique ID of the current dataset for
#' a given population and species.
#' @param DD Character specifying how to account for density dependence by
#' including population size as additional explanatory in the SEM. Possibilities
#' are: 'none' - no inclusion of density dependence, 'n_effectGR' - effect
#' of population size on growth rate, 'n_effectD' - effect of population size
#' on demographic rate, 'n_effectDGR' - effect of population size on demographic
#' rate and growth rate. Only options 'n_effectGR' and 'none' are possible for
#' simple structure of SEM, in which climate effect on GR i sonly mediated by trait
#' (i.e. \code{simpleSEM = TRUE}).
#' @param weight Logical (TRUE/FALSE) specifying whether to use inverse of
#' variances for traits and demographic rates (if \code{simpleSEM = FALSE})
#' as weights in the respective models. Defaults to FALSE.
#' @param correlation Logical (TRUE/FALSE) specifying whether to include
#' temporal autocorrelation with AR(1) structure. Defaults to FALSE.
#' @param Trait Logical (TRUE/FALSE) specifying whether to include
#' the direct impact of trait on GR in the SEM fitted using a complex
#' structure, i.e. climate effect on GR is mediated by both trait and
#' demographic rate. Is applicable only if \code{simpleSEM = FALSE}.
#' Defaults to FALSE.
#' @param simpleSEM Logical (TRUE/FALSE) specifying that the simple SEM
#' structure will be fitted in which climate effect on GR is only mediated
#' by trait, demographic rates are not considered. Defaults to FALSE.
#' @param ... Additional parameters for the function \code{\link{fit_mod}}
#'
#' @export
#' @importFrom magrittr "%>%"
#'
#' @return A list of length two, of which one is the fitted object
#'  of class psem (see \code{\link[piecewiseSEM]{psem}} for more details)
#'  and the second one is a data frame with the results on variance
#'  partitioning.
#'
#' @examples
#'
#'  # prepare the data to fit the model to one study only
#'  subs <- droplevels(dataSEM[dataSEM$ID == 7, ])
#'  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1))
#'  consec_yrs <- merge(full_NA, subs, by = 'Year', all= TRUE)
#'
#'  data_GR <- consec_yrs %>%
#'  dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
#'  dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
#'  dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
#'                !is.na(Demog_rate_mean) & !is.na(Pop_mean)) %>%
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

fit_mod <- function(biol_data, ID,
                    DD = 'none',
                    weight = FALSE,
                    correlation = FALSE,
                    Trait = FALSE,
                    simpleSEM = FALSE,
                    ...){

  dat <- droplevels(biol_data[biol_data$ID == ID, ])


  # formulas
  if(simpleSEM){
    if(DD == 'n_effectGR'){
      formGR <<- 'GR ~ det_Clim  + Pop_mean + Trait_mean'
    }
    if(DD == 'n_effectD'){
      stop("If simpleSEM = TRUE, 'DD' argument must be \"n_effectGR\" or \"none\"")
    }
    if(DD == 'n_effectDGR'){
      stop("If simpleSEM = TRUE, 'DD' argument must be \"n_effectGR\" or \"none\"")
    }
    if(DD == 'none'){
      formGR <<- 'GR ~ det_Clim + Trait_mean'
    }
  } else {
  if(Trait){

    if(DD == 'n_effectGR'){
      formGR <<- 'GR ~ det_Clim + Demog_rate_mean + Pop_mean + Trait_mean'
      formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean'
    }
    if(DD == 'n_effectD'){
      formGR <<- 'GR ~ det_Clim + Demog_rate_mean + Trait_mean'
      formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean + Pop_mean'
    }
    if(DD == 'n_effectDGR'){
      formGR <<- 'GR ~ det_Clim + Demog_rate_mean + Pop_mean + Trait_mean'
      formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean + Pop_mean'
    }
    if(DD == 'none'){
      formGR <<- 'GR ~ det_Clim + Demog_rate_mean + Trait_mean'
      formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean'
    }

  } else {
    if(DD == 'n_effectGR'){
      formGR <<- 'GR ~ det_Clim + Demog_rate_mean + Pop_mean'
      formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean'
    }
    if(DD == 'n_effectD'){
      formGR <<- 'GR ~ det_Clim + Demog_rate_mean'
      formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean + Pop_mean'
    }
    if(DD == 'n_effectDGR'){
      formGR <<- 'GR ~ det_Clim + Demog_rate_mean + Pop_mean'
      formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean + Pop_mean'
    }
    if(DD == 'none'){
      formGR <<- 'GR ~ det_Clim + Demog_rate_mean'
      formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean'
    }

  }
  }

  formTrait <<- 'Trait_mean ~ det_Clim'

  if (weight) {
    if (! simpleSEM) {
      if (sum(is.na(dat$Demog_rate_SE)) == nrow(dat)){
        dat$weights_DemRate <- 1
      } else {
        if(! any(is.na(dat$Demog_rate_SE))){
          # replace SE of 0 with min values observed, otherwise weights are Inf
          if(sum(dat$Demog_rate_SE == 0, na.rm = TRUE) != 0){
            dat$Demog_rate_SE[dat$Demog_rate_SE == 0] <-
              min(dat$Demog_rate_SE[dat$Demog_rate_SE != 0], na.rm = T)
          }
          dat$weights_DemRate <- 1 / dat$Demog_rate_SE^2
        } else {
          if(sum(dat$Demog_rate_SE == 0, na.rm = TRUE) != 0){
            dat$Demog_rate_SE[dat$Demog_rate_SE == 0] <-
              min(dat$Demog_rate_SE[dat$Demog_rate_SE != 0], na.rm = T)
          }
          dat$weights_DemRate[is.na(dat$Demog_rate_SE)] <- stats::median(dat$Demog_rate_SE, na.rm = T)  ## setting to 1 only makes sense if the analyses are run on standardized data
          dat$weights_DemRate[!is.na(dat$Demog_rate_SE)] <- 1 / dat$Demog_rate_SE[!is.na(dat$Demog_rate_SE)]^2
        }
      }
    }
    dat$weights_Trait <-  1 / dat$Trait_SE^2

  } else {
    dat$weights_DemRate <- rep(1, nrow(dat))
    dat$weights_Trait <- rep(1, nrow(dat))
  }


  if(correlation){

    dat <<- dat  ## for autocorrelation, otherwise gls does not wok within psem
    if (simpleSEM){
      partR2 <- partitionR2(data = dat)
      script <- paste0("models_list <- piecewiseSEM::psem(nlme::gls(",
                       formGR, ", correlation = nlme::corAR1(form = ~ Year | ID),
                     method = 'REML', data = dat),
                     nlme::gls(", formTrait, ",
                     weights = nlme::varFixed(~1/weights_Trait),
                     correlation = nlme::corAR1(form = ~ Year | ID),
                     method = 'REML', data = dat),
                     data = dat)")
    } else {
      script <- paste0("models_list <- piecewiseSEM::psem(nlme::gls(",
                       formGR, ", correlation = nlme::corAR1(form = ~ Year | ID),
                     method = 'REML', data = dat),
                     nlme::gls(", formDemRate, ",
                     weights = nlme::varFixed(~1/weights_DemRate),
                     correlation = nlme::corAR1(form = ~ Year | ID),
                     method = 'REML', data = dat),
                     nlme::gls(", formTrait, ",
                     weights = nlme::varFixed(~1/weights_Trait),
                     correlation = nlme::corAR1(form = ~ Year | ID),
                     method = 'REML', data = dat),
                     data = dat)")
    }


   suppressWarnings(eval(parse(text = script)))

  } else {  ## no autocorrelation

    if(simpleSEM) {
  weights_Trait <- NULL
      models_list <- piecewiseSEM::psem(
        stats::lm(stats::as.formula(formGR),
                  data = dat, ...),
        stats::lm(stats::as.formula(formTrait),
                  weights = weights_Trait,
                  data = dat, ...),
        data = dat)

    } else {
      weights_DemRate <- NULL
      models_list <- piecewiseSEM::psem(
        stats::lm(stats::as.formula(formGR),
                  data = dat, ...),
        stats::lm(stats::as.formula(formDemRate),
                  weights = weights_DemRate,
                  data = dat, ...),
        stats::lm(stats::as.formula(formTrait),
                  weights = weights_Trait,
                  data = dat, ...),
        data = dat)
    }
  }
  return(tibble::tibble(models_list = list(models_list), partR2 = list(partR2)))
}
