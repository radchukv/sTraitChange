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
#' including population size as additional explanatory in the model. Possibilities
#' are: 'none' - no inclusion of density dependence, 'n_effectGR' - effect
#' of population size on growth rate, 'n_effectD' - effect of population size
#' on demographic rate, 'n_effectDGR' - effect of population size on demographic
#' rate and growth rate.
#' @param weights Logical (TRUE/FALSE) specifying whether to use inverse of
#' variances for traits and demographic rates as weights in the respective
#' models.
#' @param correlation Logical (TRUE/FALSE) specifying whether to include
#' temporal autocorrelation with AR(1) structure.
#'
#' @export
#'
#' @return A fitted object of class psem (see
#' \code{\link[piecewiseSEM]{psem}} for more details
#'
#' @examples
#' biol_dat <- read.csv('./data-raw/Closed_Data_26_02.csv')
#' biol_stand <- convert_JulianDay(biol_data = biol_dat)
#' # example for ID = 1
#' subs <- droplevels(biol_stand[biol_stand$ID == 1, ])
#' test_fSEM <- readRDS(paste0('./output_climwin/', subs$ID[1], '_',
#'                     subs$Species[1], '_', subs$Location[1],
#'                     '_', subs$Trait[1], '_ForSEM',  '.RDS'))
#' dat <- test_fSEM$data_res[[1]]
#' ## impute the DR
#' dat <- impute_ma(data = dat, column = 'Demog_rate_mean')
#' dat <- impute_median(data = dat, column = 'Demog_rate_SE')
#' ## impute the Pop Size
#' dat <- impute_ma(data = dat, column = 'Pop_mean')
#' if(unique(subs$Count) == 'N'){
#' dat <- impute_median(data = dat, column = 'Pop_SE')
#' }
#' ## calculate GR
#' data_GR <- dat %>%
#' dplyr::mutate(., Pop_mean_lag = lag(Pop_mean)) %>%
#' dplyr::mutate(., GR = Pop_mean / Pop_mean_lag) %>%
#' dplyr::filter(., !is.na(GR))
#'
#' # fit the model
#' test <- fit_mod(biol_data = data_GR, ID =1,
#'                 DD = FALSE, weight = FALSE,
#'                 correlation = FALSE)  ## does not work for corr = TRUE yet

fit_mod <- function(biol_data, ID,
                    DD = 'none',
                    weight = FALSE,
                    correlation = FALSE,
                    ...){

  data <- droplevels(biol_data[biol_data$ID == ID, ])

  ## this will have to be taken care of at the data cleaning step...
 #  if(! any(is.na(data$Demog_rate_SE))){
 #  data <- data[data$Demog_rate_SE != 0, ]
 # }


  # formulas
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
    formGR <<- 'GR ~ det_Clim + Demog_rate_mean'   ## gls have problmes with environment...
    formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean'
  }

  formTrait <<- 'Trait_mean ~ det_Clim'

  # weights - what to do if only a few SEs are missing??? (we need weights)
  # + if all SEs are missing weights should be set to 1s
  if (weight) {
    if (sum(is.na(data$Demog_rate_SE)) == nrow(data)){
      data$weights_DemRate <- 1
    } else {
      if(! any(is.na(data$Demog_rate_SE))){
        # replace SE of 0 with min values observed, otherwise weights are Inf
        if(sum(data$Demog_rate_SE == 0) != 0){
          data$Demog_rate_SE[data$Demog_rate_SE == 0] <-
            min(data$Demog_rate_SE[data$Demog_rate_SE != 0], na.rm = T)
          }
        data$weights_DemRate <- 1 / data$Demog_rate_SE^2
      } else {
        if(sum(data$Demog_rate_SE == 0) != 0){
          data$Demog_rate_SE[data$Demog_rate_SE == 0] <-
            min(data$Demog_rate_SE[data$Demog_rate_SE != 0], na.rm = T)
        }
        data$weights_DemRate[is.na(data$Demog_rate_SE)] <- median(data$Demog_rate_SE, na.rm = T)  ## setting to 1 only makes sense if the analyses are run on standardized data
        data$weights_DemRate[!is.na(data$Demog_rate_SE)] <- 1 / data$Demog_rate_SE[!is.na(data$Demog_rate_SE)]^2
      }
    }
    data$weights_Trait <-  1 / data$Trait_SE^2  ## do same as for SE on Dem.rates?

  } else {
    data$weights_DemRate <- rep(1, nrow(data))
    data$weights_Trait <- rep(1, nrow(data))
  }

  #data <<- data  ## for autocorrelation, otherwise gls does not wok within psem
  if(correlation){
    models_list <- piecewiseSEM::psem(
      nlme::gls(stats::as.formula(formGR),  ## maybe try to go not with gls but another model that allows for weights and autocor
                data = data,
                method = 'REML', ...),
      nlme::gls(stats::as.formula(formDemRate),
                weights = nlme::varFixed(~1/weights_DemRate),
                data = data,
                method = 'REML', ...),
      nlme::gls(stats::as.formula(formTrait),
                weights = nlme::varFixed(~1/weights_Trait),
                data = data,
                method = 'REML', ...),
      data = data)
  } else {
    models_list <- piecewiseSEM::psem(
      stats::lm(stats::as.formula(formGR),
                data = data,
                method = 'qr', ...),
      stats::lm(stats::as.formula(formDemRate),
                weights = weights_DemRate,
                data = data,
                method = 'qr', ...),
      stats::lm(stats::as.formula(formTrait),
                weights = weights_Trait,
                data = data,
                method = 'qr', ...),
      data = data)
  }
  return(models_list)

}
