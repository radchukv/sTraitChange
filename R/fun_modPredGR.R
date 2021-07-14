#' Predict GR by climate and trait from the fitted linear or mixed-effects model
#'
#' \code{fun_modPredGR} predicts GR to the new data, whose columns are specified by the
#' user. The GR is predicted using the fitted linear or mixed-effects
#' model with trait and climate as predictors.
#'
#' @param mods A fitted mixed-effects or linear model with GR as a response and climate
#' and trait (as well as their quadratic effects, if desirable) as fixed effects. The linear
#' models are returned by the function \code{fit_nonLinGR}.
#' @param ME A Boolean specifying whether the supplied \code{mods} is a mixed-effects
#' model (TRUE) or a linear model (FALSE).
#' @param Year The value(s) of year to be used for predicting GR, is needed only for
#' the mixed-effects model.
#' @param det_Clim The value(s) of det_Clim to be used for predicting GR.
#' @param Trait_mean The value(s) of Trait_mean to be used for predicting GR.
#' @param Pop_mean The value(s) of Pop_mean to be used for predicting GR.
#' @inheritParams plot_relation
#'
#' @export
#'
#' @return A data frame with the requested values for all predictor variables and
#' the predicted GR values.
#'
#' ## Still add examples
#'
fun_modPredGR <- function(mods, ID,
                          Trait_categ = 'Phenological',
                          det_Clim = 0,
                          Year = 2000,
                          Trait_mean = seq(-4, 4, by = 0.25),
                          Pop_mean = 0,
                          ME = FALSE){


  data_new <- tidyr::expand_grid(det_Clim, Trait_mean, Pop_mean)
  data_new$det_Clim2 <- data_new$det_Clim ^ 2
  data_new$Trait_mean2 <- data_new$Trait_mean ^ 2
  data_new$ID <- ID


  if (! ME){
    data_new$Trait_Categ <- Trait_categ
      data_new$PredGR <- predict(mods, newdata = data_new)
  } else {
    data_new$Year <- Year
    data_new$predGR <- as.numeric(predict(mods, newdata = data_new, re.form = NA)) ## we use NA because this is already a meta-analysis
    ## add SE
    data_new$PredSE <- sqrt(spaMM::get_predVar(mods, newdata = data_new))
    data_new$lwr <- data_new$predGR + stats::qnorm(0.025) * data_new$PredSE
    data_new$upr <- data_new$predGR + stats::qnorm(0.975) * data_new$PredSE

  }

  return(data_new)
}
