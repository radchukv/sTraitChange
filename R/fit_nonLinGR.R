#' Compare fit of mixed-effects models used to predict GR
#' by non-linear effects of climate and trait
#'
#' \code{fit_nonLinGR} prepares data in a way analogous to that used to fit SEM, fits
#' four mixed-effects models differing in how they incorporate quadratic effect of
#' climate and / or trait on GR, and compares these models (using AIC)
#'
#' @inheritParams fit_SEM
#'
#' @export
#'
#' @return If all models fit well, returns a tibble, that, in addition
#' to the fields contained in the object returned by \code{fit_mod_nonL},
#' also contains the column (ModSt_minAIC) specifying the model with the
#' minimum AIC, and the column (minAICmod). The data frame, in
#' addition to the columns returned by \code{fit_mod_nonL}, also contains
#' columns ModSt_minAIC (a character specifying the 'best' model, i.e. the
#' one with the lowest AIC), BestMod2 (the second-best model as judged by AIC),
#' DeltAIC_best2 (a numeric giving delta AIC between the best and the
#' second-best model), DeltAIC_lin (a numeric specifying delta AIC of the
#' best model with the full linear model), and the coefficient estimates
#' for linear and quadratic (if applicable) effects of climate on trait,
#' for the best model.
#' If at least one model fails to fit, returns a warning.
#'
#'
#' ## Still add examples
#'
fit_nonLinGR <- function(biol_data, ID, out_nonLinM,
                         correlation = FALSE,
                         standardize = FALSE){

  # select one study
  subs <- droplevels(biol_data[biol_data$ID == ID, ])

  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1))
  consec_yrs <- merge(full_NA, subs, by = 'Year', all= T)

  message(paste('Currently fitting the model for study',
                ID, 'for species', subs$Species[1],
                'in', subs$Location[1], 'for', subs$Trait[1]))

  ## calculate GR
  data_GR <- consec_yrs %>%
    dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
    dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
    dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
                    !is.na(Demog_rate_mean) & !is.na(Pop_mean)) %>%
    dplyr::mutate(det_Clim = stats::resid(stats::lm(Clim ~ Year,
                                                    data = .)))


  if(standardize){
    data_GR <- data_GR %>%
      dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                    Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                    det_Clim = scale(det_Clim),
                    Trait_mean = scale(Trait_mean),
                    Demog_rate_mean = scale(Demog_rate_mean),
                    Pop_mean = scale(Pop_mean),
                    GR = scale(GR))
  }

  ## now call a function fitting a model (depending on the options: - autocor / no)
  fittedMod <- fit_mod_nonL(biol_data = data_GR, ID = ID,
                            correlation = correlation)

  if (! class(fittedMod)[1] == 'character'){

    ## here compare AICc and get the one with min
    minAIC_mod <- fittedMod$data[[1]]$Model[fittedMod$data[[1]]$AIC == min(fittedMod$data[[1]]$AIC)]
    fittedMod$minAICmod <- fittedMod[[paste0(minAIC_mod, '_mod')]]
    fittedMod$ModSt_minAIC <- minAIC_mod
    fittedMod$data[[1]]$ModSt_minAIC <- minAIC_mod

    ## getting the deltaAIC to the next best model and to the linear model
    nextbestMod <- fittedMod$data[[1]]$Model[fittedMod$data[[1]]$AIC == min(fittedMod$data[[1]]$AIC[fittedMod$data[[1]]$Model != fittedMod$data[[1]]$ModSt_minAIC])]
    fittedMod$data[[1]]$BestMod2 <- nextbestMod
    fittedMod$data[[1]]$DeltAIC_best2 <- min(fittedMod$data[[1]]$AIC[fittedMod$data[[1]]$Model != fittedMod$data[[1]]$ModSt_minAIC]) - min(fittedMod$data[[1]]$AIC)
    fittedMod$data[[1]]$DeltAIC_lin <- fittedMod$data[[1]]$AIC[fittedMod$data[[1]]$Model == 'fullLin'] - min(fittedMod$data[[1]]$AIC)


    ## also add the coefficient estimates to the data frame within the tibble
    fittedMod$data[[1]]$det_Clim <- coef(fittedMod$minAICmod[[1]])['det_Clim']
    fittedMod$data[[1]]$Trait_mean <- coef(fittedMod$minAICmod[[1]])['Trait_mean']
    fittedMod$data[[1]]$Pop_mean <- coef(fittedMod$minAICmod[[1]])['Pop_mean']
    fittedMod$data[[1]]$det_Clim2 <- coef(fittedMod$minAICmod[[1]])['det_Clim2']
    fittedMod$data[[1]]$Trait_mean2 <- coef(fittedMod$minAICmod[[1]])['Trait_mean2']

    # perhaps it may be useful to write output per model to the disk, for now left out

    return(fittedMod)
  } else { ## i.e. the model did not fit with fit_mod_nonL for whatever reason
    warning(cat('error when fitting the model \n'))
  }
}
