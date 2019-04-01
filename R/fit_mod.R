#' Fit structural equation model with requested specifications
#'
#' \code{fit_mod} fits structural equation model with
#' requested model structure for the specified study ID
#'
#' @param biol_data Data frame with climate, trait and demographic rate
#' data for a given population and species.
#' @param ID Numeric giving a unique ID of the current dataset for
#' a given population and species.
#' @param DD Logical (TRUE/FALSE) specifying whether to account for density
#' dependence by including population size as additional explanatory of
#' population growth rate.
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
#'                 DD = FALSE, weights = FALSE,
#'                 correlation = FALSE)  ## does not work for corr = TRUE yet

fit_mod <- function(biol_data, ID,
                    DD = FALSE,
                    weights = FALSE,
                    correlation = FALSE,
                    ...){

  data <- droplevels(biol_data[biol_data$ID == ID, ])
  # formulas
  if(DD){
    formGR <<- 'GR ~ Clim + Demog_rate_mean + Pop_mean'
  }else{
    formGR <<- 'GR ~ Clim + Demog_rate_mean'   ## gls have problmes with environment...
  }                   ## should probably contact the developers...

  formDemRate <<- 'Demog_rate_mean ~ Clim + Trait_mean'
  formTrait <<- 'Trait_mean ~ Clim'

  # weights
  if (weights) {
    data$weights_DemRate <- 1 / data$Demog_rate_SE^2
    data$weights_Trait <-  1 / data$Trait_SE^2
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
