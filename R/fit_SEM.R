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
#' # for this example we will be saving the output plots in the tempdir folder
#' # let's check its path and content before we run the example
#' message('Temporary directory is located at', tempdir())
#' message('Contents of the temporary directory', list.files(tempdir()))
#' mod_SEM <- fit_SEM(biol_data = dataSEM, ID = 1,
#'                 out_SEM = tempdir(), # attention: for this example we write the data to
#'                                      # a temporary directory, to check its location type tempdir()
#'                 DD = 'n_effectGR', weight = TRUE,
#'                 correlation = TRUE,
#'                 standardize = TRUE,
#'                 Trait = FALSE,
#'                 simpleSEM = TRUE)
#' message('Contents of the temporary directory after running fit_SEM()',
#' list.files(tempdir()))
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
    dplyr::mutate(Pop_mean_lag = c(.data$Pop_mean[-1], NA)) %>%
    dplyr::mutate(GR = log(.data$Pop_mean_lag / .data$Pop_mean)) %>%
    dplyr::filter(!is.na(.data$GR) & !is.na(.data$Trait_mean) &
                    !is.na(.data$Demog_rate_mean) & !is.na(.data$Pop_mean)) %>%
    dplyr::mutate(det_Clim = stats::resid(stats::lm(Clim ~ Year,
                                                    data = .))) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.array), as.numeric))


  # exploratory plots
  #message(paste(out_SEM))
  path <- paste0(out_SEM, '/', data_GR$ID[1], '_',
                 data_GR$Species[1], '_', data_GR$Location[1],
                 '_', data_GR$Trait[1], '_relations.pdf')
  # message(paste(path))

  if (requireNamespace("psych", quietly = TRUE)) {
  grDevices::pdf(path)
  dat_sub <- data_GR %>%
    dplyr::select(.data$Clim, .data$det_Clim, .data$Year,
                  .data$Trait_mean, .data$Demog_rate_mean,
                  .data$Pop_mean, .data$GR)
  psych::pairs.panels(dat_sub, ellipses = FALSE, hist.col = 'grey', lm = TRUE)
  graphics::mtext(paste0('Demographic rate is ', unique(dat_sub$Demog_rate)), side = 3,
        line = 3)
  graphics::mtext(paste0('Trait is ', unique(dat_sub$Trait_Categ_det)), side = 1,
        line = 4)
  grDevices::dev.off()
  } else {
    message("to be able to produce a pdf with the raw data, you must first run install.packages('psych')!")
  }


  if(standardize){
    data_GR <- data_GR %>%
      dplyr::mutate(Trait_SE = .data$Trait_SE / stats::sd(.data$Trait_mean, na.rm = T),
                    Demog_rate_SE = .data$Demog_rate_SE / stats::sd(.data$Demog_rate_mean, na.rm = T),
                    det_Clim = scale(.data$det_Clim),
                    Trait_mean = scale(.data$Trait_mean),
                    Demog_rate_mean = scale(.data$Demog_rate_mean),
                    Pop_mean = scale(.data$Pop_mean),
                    GR = scale(.data$GR)) %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.array), as.numeric))



    #message(paste(out_SEM))
    path <- paste0(out_SEM, '/', data_GR$ID[1], '_',
                   data_GR$Species[1], '_', data_GR$Location[1],
                   '_', data_GR$Trait[1], '_z_score_relations.pdf')
    #message(paste(path))

    if (requireNamespace("psych", quietly = TRUE)) {
    grDevices::pdf(path)
    dat_sub <- data_GR %>%
      dplyr::select(.data$Clim, .data$det_Clim, .data$Year,
                    .data$Trait_mean, .data$Demog_rate_mean,
                    .data$Pop_mean, .data$GR)
    psych::pairs.panels(dat_sub, ellipses = FALSE, hist.col = 'grey', lm = TRUE)
    graphics::mtext(paste0('Demographic rate is ', unique(dat_sub$Demog_rate)), side = 3,
          line = 3)
    graphics::mtext(paste0('Trait is ', unique(dat_sub$Trait_Categ_det)), side = 1,
          line = 4)
    grDevices::dev.off()
    } else {
      message("to be able to produce a pdf with the raw data, you must first run install.packages('psych')!")
    }

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
  SEM_results <- get_res_SEM(SEM_mod$models_list[[1]])

  SEM_results$DD <- DD
  SEM_results$weights <- weight
  SEM_results$corr <- correlation
  SEM_results$Tr <- Trait

  # write output
  path <- paste0(out_SEM, '/', SEM_results$ID[1], '_',
                 SEM_results$Species[1], '_', SEM_results$Location[1],
                 '_', SEM_results$Trait[1], '_',
                 SEM_results$Demog_rate[1], '_ResultsSEM',  '.RDS')
 # message(paste(path))
  saveRDS(object = SEM_results,
          file = path)
  out <- tibble::tibble(SEM_results = list(SEM_results),
                       partR2 = list(SEM_mod$partR2[[1]]))
  return(out)
}
