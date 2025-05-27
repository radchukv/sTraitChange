#' Partition R2 in the full model
#'
#' \code{partitionR2} partitions the variation in the data explained by
#' fixed effects of the fitted GR-model (the model to explain population
#' growth rate). Currently only works for the model that includes
#' autocorrelation in residual structure and simpleSEM  option set to TRUE
#'
#' @param data Data frame for a single study that contains all the
#' covariates needed to fit a single model
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export
#'
#' @return A data frame with four columns, each for a variation explained by:
#' 1) population size, 2)climate, 3) trait, 4) the full model.
#'
#' @examples
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
#'  dat <- droplevels(data_GR[data_GR$ID == 7, ])
#'
#' # partition the variation
#' partR2_study7 <- partitionR2(data = dat)
#'
partitionR2 <- function(data){
# fit full model
mod_full <- nlme::gls(GR ~ det_Clim + Trait_mean + Pop_mean,
                      correlation = nlme::corAR1(form = ~ Year),
          method = 'ML', data = data)
# fit each submodel without the one fixed-effect covariate in turn
mod_noPop <- nlme::gls(GR ~ det_Clim + Trait_mean,
                       correlation = nlme::corAR1(form = ~ Year),
                       method = 'ML', data = data)
mod_null <- nlme::gls(GR ~ 1,
                      method = 'ML', data = data)
mod_noClim <- nlme::gls(GR ~ Pop_mean + Trait_mean,
                        correlation = nlme::corAR1(form = ~ Year),
                        method = 'ML', data = data)
mod_noTrait <- nlme::gls(GR ~ Pop_mean + det_Clim,
                        correlation = nlme::corAR1(form = ~ Year),
                        method = 'ML', data = data)
PopsizeR2_lik <- rr2::R2_lik(mod_full, mod_noPop)
ClimateR2_lik <- rr2::R2_lik(mod_full, mod_noClim)
TraitR2_lik <- rr2::R2_lik(mod_full, mod_noTrait)
FullR2_lik <- rr2::R2_lik(mod_full, mod_null)
dat_r2 <- data.frame(PopsizeR2_lik = PopsizeR2_lik,
                     ClimateR2_lik = ClimateR2_lik,
                     TraitR2_lik = TraitR2_lik,
                     FullR2_lik = FullR2_lik)
return(dat_r2)
}
