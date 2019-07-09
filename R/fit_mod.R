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

  dat <- droplevels(biol_data[biol_data$ID == ID, ])

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
    formGR <<- 'GR ~ det_Clim + Demog_rate_mean'
    formDemRate <<- 'Demog_rate_mean ~ det_Clim + Trait_mean'
  }

  formTrait <<- 'Trait_mean ~ det_Clim'


  if (weight) {
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
        dat$weights_DemRate[is.na(dat$Demog_rate_SE)] <- median(dat$Demog_rate_SE, na.rm = T)  ## setting to 1 only makes sense if the analyses are run on standardized data
        dat$weights_DemRate[!is.na(dat$Demog_rate_SE)] <- 1 / dat$Demog_rate_SE[!is.na(dat$Demog_rate_SE)]^2
      }
    }
    dat$weights_Trait <-  1 / dat$Trait_SE^2  ## do same as for SE on Dem.rates?

  } else {
    dat$weights_DemRate <- rep(1, nrow(dat))
    dat$weights_Trait <- rep(1, nrow(dat))
  }


  if(correlation){

<<<<<<< HEAD
    dat <<- dat  ## for autocorrelation, otherwise gls does not wok within psem
=======
    # dat <<- dat  ## for autocorrelation, otherwise gls does not wok within psem
>>>>>>> cd735a441f1d66bc99af956c15596b7ed3d6d7fd
    # formGR <<- formGR
    # formDemRate <<- formDemRate
    # formTrait <<- formTrait


    # trying with Alex's suggestion
    script <- paste0("models_list <- piecewiseSEM::psem(nlme::gls(",
                     formGR, ", method = 'REML', data = dat),
                     nlme::gls(", formDemRate, ",
                     weights = nlme::varFixed(~1/weights_DemRate),
                     method = 'REML', data = dat),
                     nlme::gls(", formTrait, ",
                     weights = nlme::varFixed(~1/weights_Trait),
                     method = 'REML', data = dat),
                     data = dat)")

    eval(parse(text = script))


    # models_list <- piecewiseSEM::psem(
    #   nlme::gls(stats::as.formula(formGR),  ## maybe try to go not with gls but another model that allows for weights and autocor
    #             data = dat,
    #             method = 'REML', ...),
    #   nlme::gls(stats::as.formula(formDemRate),
    #             weights = nlme::varFixed(~1/weights_DemRate),
    #             data = dat,
    #             method = 'REML', ...),
    #   nlme::gls(stats::as.formula(formTrait),
    #             weights = nlme::varFixed(~1/weights_Trait),
    #             data = dat,
    #             method = 'REML', ...),
    #   data = dat)
  } else {
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
  return(models_list)

}
