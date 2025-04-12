#' Partition R2 in the full model
#'
#' \code{partitionR2} partitions the variation in the data explained by
#' the full model due by its fixed-effects
#'
#' @param data Data frame for a single study that contains all the
#' covariates needed to fit a single model
#' @inheritParams fit_mod
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export
#'
#' @return A data frame with four columns, each for a variation explained by:
#' 1) population size, 2)climate, 3) trait, 4) the full model.
#'
#' @examples
#' # still to do
#'
# think - I probably have to build several if-parts depending on whether
# the autocor is included or not... BUT as I am currently fittign one specific
# mdoel, maybe it does not really matter???
partitionR2 <- function(data = dat){
# fit full model
mod_full <- nlme::gls(GR ~ det_Clim + Trait_mean + Pop_mean,
                      correlation = nlme::corAR1(form = ~ Year | ID),
          method = 'REML', data = dat)
# fit each submodel without theone fixed-effect covariate in turn
mod_noPop <- nlme::gls(GR ~ det_Clim + Trait_mean,
                       correlation = nlme::corAR1(form = ~ Year | ID),
                       method = 'REML', data = dat)
mod_null <- nlme::gls(GR ~ 1,
                      method = 'REML', data = dat)
mod_noClim <- nlme::gls(GR ~ Pop_mean + Trait_mean,
                        correlation = nlme::corAR1(form = ~ Year | ID),
                        method = 'REML', data = dat)
mod_noTrait <- nlme::gls(GR ~ Pop_mean + det_Clim,
                        correlation = nlme::corAR1(form = ~ Year | ID),
                        method = 'REML', data = dat)
PopsizeR2 <- R2_lik(mod_full, mod_noPop)
ClimateR2 <- R2_lik(mod_full, mod_noClim)
TraitR2 <- R2_lik(mod_full, mod_noTrait)
FullR2 <- R2_lik(mod_full, mod_null)  # 0.1 explained in total, of which 0.45
dat_r2 <- data.frame(PopsizeR2 = PopsizeR2, ClimateR2 = ClimateR2,
                     TraitR2 = TraitR2, FullR2 = FullR2)
return(dat_r2)
}
