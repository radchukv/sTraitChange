#' Test whether nonlinear effect fits better across studies
#'
#' \code{test_nonlin} fits the model with and without quadratic fixed effect
#' across the studies. The function fits the model using REML with and
#' without the AR1 auto-regressive term. It then compares the AIC of these
#' two models. The model with best selected random structure is then
#' refitted using ML and compared to a null model.
#' In a next stepm, in case fixed quadratic effect is supported by the data,
#' this "random intercept" model is compared with the random slope (where
#' random slope is just included in a random effect as a linear effect),
#' using REML.
#' All models are fitted using the package \pkg{\link[spaMM]{spaMM}}.
#'
#' @inheritParams parameters_definition
#'
#' @return A list containing the best full fitted model and the result of
#'   \code{\link[spaMM]{anova.HLfit}} comparing this model and a null model.
#' @export
#'
#' ## still add the examples
#'
test_nonlin <- function(data,
                        formula_full  = '~ 1',
                        formula_null  = '~ 1',
                        rand_sl = '+ (det_Clim|ID)',
                        ...) {

  mod_full_fix <- spaMM::fitme(stats::as.formula(formula_full),
                               data = data,
                               method = "REML", ...)

  formula_full_AR <- paste(formula_full, '+ AR1(1|Year)')

  mod_full_AR <- spaMM::fitme(stats::as.formula(formula_full_AR),
                              data = data,
                              method = "REML", ...)

  fix_best <- spaMM::AIC.HLfit(mod_full_fix)[[1]] < spaMM::AIC.HLfit(mod_full_AR)[[1]]

  out <- list()

  if (!fix_best) {
    formula_full <- paste(formula_full, '+ AR1(1|Year)')
    formula_null <- paste(formula_null, '+ AR1(1|Year)')
    out$model <- mod_full_AR
  } else {
    out$model <- mod_full_fix
  }

  mod_full <- spaMM::fitme(stats::as.formula(formula_full),
                           data = data,
                           method = "ML", ...)

  mod_null <- spaMM::fitme(stats::as.formula(formula_null),
                           data = data,
                           method = "ML", ...)

  out$anova <- spaMM::anova.HLfit(mod_full, mod_null)


  if(out$anova[[3]]$p_value < 0.05){
    formula_randInterc <- formula_full
    formula_randSlope <- paste(formula_full, rand_sl)
    mod_randInterc <- spaMM::fitme(stats::as.formula(formula_randInterc),
                                   data = data,
                                   method = "REML", ...)

    mod_randSlope <- spaMM::fitme(stats::as.formula(formula_randSlope),
                                  data = data,
                                  method = "REML", ...)
    out$AIC_randInterc <- spaMM::AIC.HLfit(mod_randInterc)[[1]]
    out$AIC_randSlope <- spaMM::AIC.HLfit(mod_randSlope)[[1]]
  }
  return(out)
}
