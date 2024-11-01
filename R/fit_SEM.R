#' Fit structural equation model and return the results
#'
#' \code{fit_SEM} fits structural equation model with
#' requested model structure for the specified study ID and returns
#' the results to be used in meta-analysis
#'
#' @param out_SEM Character specifying the library on the path where
#' the results of structural equation analysis are stored.
#' @param standardize A boolean specifying whether the variables have to be
#'  standardized (z-scored) so as to lead to standardized path coefficients.
#'  Defaults to FALSE.
#' @inheritParams fit_mod
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export
#'
#' @return A tibble with 33 columns that represent the metadata for a
#' given study ID and the results of SEM: Cstat, dTable, R2 of all
#' equations and path coefficients, as well as requested settings
#' for fitting SEM (whether weights, autocorrelation and traits are
#' included in the model).
#'
#' @examples
#' # fit SEM
#' mod_SEM <- fit_SEM(biol_data = dataSEM, ID = 1,
#'                 out_SEM = 'output_SEM_all',
#'                 DD = 'n_effectGR', weight = TRUE,
#'                 correlation = TRUE,
#'                 standardize = TRUE,
#'                 Trait = FALSE,
#'                 simpleSEM = TRUE)
#'
fit_SEM <- function(biol_data, ID, out_SEM,
                    DD = 'none',
                    weight = FALSE,
                    correlation = FALSE,
                    standardize = FALSE,
                    Trait = FALSE,
                    simpleSEM = FALSE){


  # select one study
  subs <- droplevels(biol_data[biol_data$ID == ID, ])

  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1))
  consec_yrs <- merge(full_NA, subs, by = 'Year', all= T)

  message(paste('Currently fitting SEM for study',
                ID, 'for species', subs$Species[1],
                'in', subs$Location[1], 'for', subs$Trait[1]))

  ## calculate GR
  data_GR <- consec_yrs %>%
    dplyr::mutate(.data, Pop_mean_lag = c(.data$Pop_mean[-1], NA)) %>%
    dplyr::mutate(.data, GR = log(.data$Pop_mean_lag / .data$Pop_mean)) %>%
    dplyr::filter(.data, !is.na(.data$GR) & !is.na(.data$Trait_mean) &
                    !is.na(.data$Demog_rate_mean) & !is.na(.data$Pop_mean)) %>%
    dplyr::mutate(det_Clim = stats::resid(stats::lm(Clim ~ Year,
                                                    data = .))) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.array), as.numeric))


  # exploratory plots
  grDevices::pdf(paste0('./', out_SEM, '/', data_GR$ID[1], '_',
             data_GR$Species[1], '_', data_GR$Location[1],
             '_', data_GR$Trait[1], '_relations.pdf'))
  psych::pairs.panels(subset(data_GR, select = c(Clim, det_Clim, Year,
                                                 Trait_mean, Demog_rate_mean,
                                                 Pop_mean, GR)),
                      ellipses = FALSE, hist.col = 'grey', lm = TRUE)
  graphics::mtext(paste0('Demographic rate is ', unique(data_GR$Demog_rate)), side = 3,
        line = 3)
  graphics::mtext(paste0('Trait is ', unique(data_GR$Trait_Categ_det)), side = 1,
        line = 4)
  grDevices::dev.off()



  if(standardize){
    data_GR <- data_GR %>%
      dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                    Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                    det_Clim = scale(det_Clim),
                    Trait_mean = scale(Trait_mean),
                    Demog_rate_mean = scale(Demog_rate_mean),
                    Pop_mean = scale(Pop_mean),
                    GR = scale(GR)) %>%
      dplyr::mutate(across(where(is.array), as.numeric))

    grDevices::pdf(paste0('./', out_SEM, '/', data_GR$ID[1], '_',
               data_GR$Species[1], '_', data_GR$Location[1],
               '_', data_GR$Trait[1], '_z_score_relations.pdf'))
    psych::pairs.panels(subset(data_GR, select = c(Clim, det_Clim, Year,
                                                   Trait_mean, Demog_rate_mean,
                                                   Pop_mean, GR)),
                        ellipses = FALSE, hist.col = 'grey', lm = TRUE)
    graphics::mtext(paste0('Demographic rate is ', unique(data_GR$Demog_rate)), side = 3,
          line = 3)
    graphics::mtext(paste0('Trait is ', unique(data_GR$Trait_Categ_det)), side = 1,
          line = 4)
    grDevices::dev.off()
  }

  ## now call a function fitting a model (depending on the options:
  ## - autocor / no, DD/no/on what element, weights /no)

  SEM_mod <- fit_mod(biol_data = data_GR, ID = ID,
                     DD = DD,
                     weight = weight,
                     correlation = correlation,
                     Trait = Trait,
                     simpleSEM = simpleSEM)


  ## then get output from the results of fitting the function
  SEM_results <- get_res_SEM(SEM_mod)

  SEM_results$DD <- DD
  SEM_results$weights <- weight
  SEM_results$corr <- correlation
  SEM_results$Tr <- Trait

  # write output
  saveRDS(object = SEM_results,
          file = paste0('./', out_SEM, '/', SEM_results$ID[1], '_',
                        SEM_results$Species[1], '_', SEM_results$Location[1],
                        '_', SEM_results$Trait[1], '_',
                        SEM_results$Demog_rate[1], '_ResultsSEM',  '.RDS'))
  return(SEM_results)
}
