#' Fit structural equation model and return the results
#'
#' \code{fit_SEM} fits structural equation model with
#' requested model structure for the specified study ID and returns
#' the results to be used in meta-analysis
#'
#' @param biol_data Data frame with climate, trait and demographic rate
#' data for a given population and species.
#' @param ID Numeric giving a unique ID of the current dataset for
#' a given population and species.
#' @param out_SEM Character specifying the library on the path where
#' the results of structural equation analysis are stored.
#' @param DD Logical (TRUE/FALSE) specifying whether to account for density
#' dependence by including population size as additional explanatory of
#' population growth rate.
#' @param weights Logical (TRUE/FALSE) specifying whether to use inverse of
#' variances for traits and demographic rates as weights in the respective
#' models.
#' @param correlation Logical (TRUE/FALSE) specifying whether to include
#' temporal autocorrelation with AR(1) structure.
#' @param standardize Logical (TRUE/FALSE) specifying whether to
#' standardize (using z scores) the variables prior to fitting the model.
#'
#' @export
#'
#' @return A tibble with 15 columns that represent the metadata for a
#' given study ID and the results of SEM: Cstat, dTable, R2 of all
#' equations and path coefficients.
#'
#' @examples
#' biol_dat <- read.csv('./data-raw/Closed_Data_26_02.csv')
#' biol_stand <- convert_JulianDay(biol_data = biol_dat)
#' # example for ID = 1
#' subs <- droplevels(biol_stand[biol_stand$ID == 25, ])
#' test_fSEM <- readRDS(paste0('./output_forSEM_temp/', subs$ID[1], '_',
#'                     subs$Species[1], '_', subs$Location[1],
#'                     '_', subs$Trait[1], '_ForSEM',  '.RDS'))
#' dat <- test_fSEM$data_res[[1]]
#' # fit the model
#' test <- fit_SEM(biol_data = dat, ID = 1,
#'                 out_SEM = 'output_SEM',
#'                 DD = FALSE, weights = FALSE,
#'                 correlation = FALSE,
#'                 standardize = FALSE)  ## does not work for corr = TRUE yet
#'
fit_SEM <- function(biol_data, ID, out_SEM,
                    DD = 'none',
                    weight = FALSE,
                    correlation = FALSE,
                    standardize = FALSE,
                    detrend = FALSE,
                    Trait = TRUE){


  # select one study
  subs <- droplevels(biol_data[biol_data$ID == ID, ])

  full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1))
  consec_yrs <- merge(full_NA, subs, by = 'Year', all= T)

  message(paste('Currently fitting SEM for study',
                ID, 'for species', subs$Species[1],
                'in', subs$Location[1], 'for', subs$Trait[1]))

  ## no imputation to be used, drop the missing values
  # ## impute the DR
  # subs <- impute_ma(data = subs, column = 'Demog_rate_mean')
  # subs <- impute_median(data = subs, column = 'Demog_rate_SE')
  # ## impute the Pop Size
  # subs <- impute_ma(data = subs, column = 'Pop_mean')
  # if(unique(subs$Count) == 'N'){
  #   subs <- impute_median(data = subs, column = 'Pop_SE')
  # }

  ## calculate GR
  data_GR <- consec_yrs %>%
    dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
    dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
    dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
                    !is.na(Demog_rate_mean) & !is.na(Pop_mean)) %>%
    dplyr::mutate(det_Clim = stats::resid(stats::lm(Clim ~ Year,
                                                    data = .)))  ## detrend clim (by year as a trend)
  ## I may need to revise this - for studies with all SE missing...
  ## so for now exclude the cases where SEs are missing


  # exploratory plots
  pdf(paste0('./', out_SEM, '/', data_GR$ID[1], '_',
             data_GR$Species[1], '_', data_GR$Location[1],
             '_', data_GR$Trait[1], '_relations.pdf'))
  psych::pairs.panels(subset(data_GR, select = c(Clim, det_Clim, Year,
                                                 Trait_mean, Demog_rate_mean,
                                                 Pop_mean, GR)),
                      ellipses = FALSE, hist.col = 'grey', lm = TRUE)
  mtext(paste0('Demographic rate is ', unique(data_GR$Demog_rate)), side = 3,
        line = 3)
  mtext(paste0('Trait is ', unique(data_GR$Trait_Categ_det)), side = 1,
        line = 4)
  dev.off()



  if(standardize){
    data_GR <- data_GR %>%
      dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = T),
                    Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = T),
                    det_Clim = scale(det_Clim),
                    Trait_mean = scale(Trait_mean),
                    Demog_rate_mean = scale(Demog_rate_mean),
                    Pop_mean = scale(Pop_mean),
                    GR = scale(GR))

    pdf(paste0('./', out_SEM, '/', data_GR$ID[1], '_',
               data_GR$Species[1], '_', data_GR$Location[1],
               '_', data_GR$Trait[1], '_z_score_relations.pdf'))
    psych::pairs.panels(subset(data_GR, select = c(Clim, det_Clim, Year,
                                                   Trait_mean, Demog_rate_mean,
                                                   Pop_mean, GR)),
                        ellipses = FALSE, hist.col = 'grey', lm = TRUE)
    mtext(paste0('Demographic rate is ', unique(data_GR$Demog_rate)), side = 3,
          line = 3)
    mtext(paste0('Trait is ', unique(data_GR$Trait_Categ_det)), side = 1,
          line = 4)
    dev.off()
  }

  if(detrend){
    data_GR$weights_Trait <-  1 / data_GR$Trait_SE^2
    data_GR <- data_GR %>%
      dplyr::mutate(Clim = scale(Clim),
                    Trait_mean = stats::resid(stats::lm(Trait_mean ~ Pop_mean,
                                                        weights = weights_Trait, data = .)) /
                      sd(stats::resid(stats::lm(Trait_mean ~ Pop_mean,
                                                weights = weights_Trait,
                                                data = .)), na.rm = TRUE),
                    Demog_rate_mean =
                      stats::resid(stats::lm(Demog_rate_mean ~ Pop_mean,
                                             data = .)) /
                      sd(stats::resid(stats::lm(Demog_rate_mean ~ Pop_mean,
                                                data = .)), na.rm = TRUE),
                    GR = stats::resid(stats::lm(GR ~ Pop_mean, data = .)) /
                      sd(stats::resid(stats::lm(GR ~ Pop_mean, data = .)), na.rm = T))
  }

  ## now call a function fitting a model (depending on the options:
  ## - autocor / no, DD/no/on what element, weights /no)

  SEM_mod <- fit_mod(biol_data = data_GR, ID = ID,
                     DD = DD,
                     weight = weight,
                     correlation = correlation,
                     Trait = Trait) ## there is another problem with gls
  ## it does not return the object if called from within psem()


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
