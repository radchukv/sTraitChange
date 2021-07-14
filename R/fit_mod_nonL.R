#' Fit mixed-effects models to predict GR by climate and trait
#'
#' \code{fit_mod_nonL} fits mixed-effects model to predict GR by
#' climate and trait, using varying structures of non-linear effects
#' for each of the predictors
#'
#' @inheritParams fit_mod
#'
#' @export
#'
#' @return If all models fit well, returns a tibble containing a
#' data frame, the study ID, the trait category, and all four
#' fitted models (full quadratic, where a quadratic effect is
#' included for both climate and trait; quadratic climate, where
#' a quadratic effect is only included for climate; quadratic traits,
#' where a quadratic effect is only included for trait; and full
#' linear model where both climate and trait have a linear effect
#' on GR). The data frame contains the study ID, the trait category,
#' and AIC and df for each of the four models.
#' If at least one model fails to fit, returns a warning.
#'
#' ## Still add examples
#'
fit_mod_nonL <- function(biol_data, ID,
                         correlation = FALSE){
  dat <- droplevels(biol_data[biol_data$ID == ID, ])
  dat$Trait_mean2 <- dat$Trait_mean ^2
  dat$det_Clim2 <- dat$det_Clim ^2


  if(correlation){

    er_fullQuad_REML <-  tryCatch(mod_fullQuad_REML <- nlme::gls(GR ~ det_Clim + det_Clim2 + Pop_mean +
                                                                   Trait_mean + Trait_mean2,
                                                                 correlation = nlme::corAR1(form = ~ Year | ID),
                                                                 method = 'REML', data = dat),
                                  error=function(e) e)
    er_fullLin_REML <- tryCatch(mod_fullLin_REML <- nlme::gls(GR ~ det_Clim  + Pop_mean +
                                                                Trait_mean,
                                                              correlation = nlme::corAR1(form = ~ Year | ID),
                                                              method = 'REML', data = dat),
                                error=function(e) e)
    er_QuadClim_REML <- tryCatch(mod_QuadClim_REML <- nlme::gls(GR ~ det_Clim + det_Clim2 + Pop_mean +
                                                                  Trait_mean,
                                                                correlation = nlme::corAR1(form = ~ Year | ID),
                                                                method = 'REML', data = dat),
                                 error=function(e) e)
    er_QuadTrait_REML <- tryCatch(mod_QuadTrait_REML <- nlme::gls(GR ~ det_Clim  + Pop_mean + Trait_mean +
                                                                    Trait_mean2,
                                                                  correlation = nlme::corAR1(form = ~ Year | ID),
                                                                  method = 'REML', data = dat),
                                  error=function(e) e)

    ## an dnow with ML, for AIC calculation
    er_fullQuad_ML <- tryCatch(mod_fullQuad_ML <- nlme::gls(GR ~ det_Clim + det_Clim2 + Pop_mean +
                                                              Trait_mean + Trait_mean2,
                                                            correlation = nlme::corAR1(form = ~ Year | ID),
                                                            method = 'ML', data = dat),
                               error=function(e) e)
    er_fullLin_ML <- tryCatch(mod_fullLin_ML <- nlme::gls(GR ~ det_Clim  + Pop_mean +
                                                            Trait_mean,
                                                          correlation = nlme::corAR1(form = ~ Year | ID),
                                                          method = 'ML', data = dat),
                              error=function(e) e)
    er_QuadClim_ML <- tryCatch(mod_QuadClim_ML <- nlme::gls(GR ~ det_Clim + det_Clim2 + Pop_mean +
                                                              Trait_mean,
                                                            correlation = nlme::corAR1(form = ~ Year | ID),
                                                            method = 'ML', data = dat),
                               error=function(e) e)
    er_QuadTrait_ML <- tryCatch(mod_QuadTrait_ML <- nlme::gls(GR ~ det_Clim  + Pop_mean + Trait_mean +
                                                                Trait_mean2,
                                                              correlation = nlme::corAR1(form = ~ Year | ID),
                                                              method = 'ML', data = dat),
                                error=function(e) e)
  } else {
    ## this is not written yet, as we anyways likely will want the model with autocor TO DO
  }

  if(! is(er_fullQuad_REML,'error') & ! is(er_fullLin_REML, 'error') &
     ! is(er_QuadClim_REML, 'error') & !is(er_QuadTrait_REML, 'error') &
     ! is(er_fullQuad_ML, 'error') & !is(er_fullLin_ML, 'error') &
     ! is(er_QuadClim_ML, 'error') & !is(er_QuadTrait_ML, 'error')){
    out_dat <- data.frame(AIC = AIC(mod_fullQuad_ML, mod_QuadClim_ML,
                                    mod_QuadTrait_ML, mod_fullLin_ML)$AIC,
                          df = AIC(mod_fullQuad_ML, mod_QuadClim_ML,
                                   mod_QuadTrait_ML, mod_fullLin_ML)$df,
                          Model = c('fullQuad', 'QuadClim', 'QuadTrait',
                                    'fullLin'),
                          ID = rep(ID, 4),
                          Trait_Categ = unique(dat$Trait_Categ))

    out_tib <- tibble::tibble(data = list(out_dat),
                              fullQuad_mod = list(mod_fullQuad_REML),
                              QuadClim_mod = list(mod_QuadClim_REML),
                              QuadTrait_mod = list(mod_QuadTrait_REML),
                              fullLin_mod = list(mod_fullLin_REML),
                              ID = ID, Trait_Categ = unique(dat$Trait_Categ))
    return(out_tib)
  } else {
    warning(cat('error when fitting the model \n'))
  }
}
