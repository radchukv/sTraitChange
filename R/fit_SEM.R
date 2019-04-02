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
#' subs <- droplevels(biol_stand[biol_stand$ID == 1, ])
#' test_fSEM <- readRDS(paste0('./output_climwin/', subs$ID[1], '_',
#'                     subs$Species[1], '_', subs$Location[1],
#'                     '_', subs$Trait[1], '_ForSEM',  '.RDS'))
#' dat <- test_fSEM$data_res[[1]]
#' # fit the model
#' test <- fit_SEM(biol_data = dat, ID = 1,
#'                 out_SEM = 'output_SEM',
#'                 DD = FALSE, weights = FALSE,
#'                 correlation = FALSE,
#'                 standardize = FALSE)  ## does not work for corr = TRUE yet

fit_SEM <- function(biol_data, ID, out_SEM,
                    DD = FALSE,
                    weights = FALSE,
                    correlation = FALSE,
                    standardize = FALSE){
  # select one study
  subs <- droplevels(biol_data[biol_data$ID == ID, ])

  ## impute the DR
  subs <- impute_ma(data = subs, column = 'Demog_rate_mean')
  subs <- impute_median(data = subs, column = 'Demog_rate_SE')
  ## impute the Pop Size
  subs <- impute_ma(data = subs, column = 'Pop_mean')
  if(unique(subs$Count) == 'N'){
    subs <- impute_median(data = subs, column = 'Pop_SE')
  }

  ## calculate GR
  data_GR <- subs %>%
    dplyr::mutate(., Pop_mean_lag = lag(Pop_mean)) %>%
    dplyr::mutate(., GR = Pop_mean / Pop_mean_lag) %>%
    dplyr::filter(., !is.na(GR))

  if(standardize){
    data_GR <- data_GR %>%
      dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean),
                    Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean),
                    Clim = scale(Clim),
                    Trait_mean = scale(Trait_mean),
                    Demog_rate_mean = scale(Demog_rate_mean),
                    Pop_mean = scale(Pop_mean),
                    GR = scale(GR))

    pdf(paste0('./', out_SEM, '/', data_GR$ID[1], '_',
               data_GR$Species[1], '_', data_GR$Location[1],
               '_', data_GR$Trait[1], '_z_score_relations.pdf'))
    psych::pairs.panels(subset(data_GR, select = c(Clim, Trait_mean,
                                                   Demog_rate_mean, Pop_mean, GR)),
                        ellipses = FALSE, hist.col = 'grey', lm = TRUE)
    mtext(paste0('Demographic rate is ', unique(data_GR$Demog_rate)), side = 3,
          line = 3)
    mtext(paste0('Trait is ', unique(data_GR$Trait_Categ_det)), side = 1,
          line = 4)
    dev.off()
     }
  # exploratory plots
  pdf(paste0('./', out_SEM, '/', data_GR$ID[1], '_',
             data_GR$Species[1], '_', data_GR$Location[1],
             '_', data_GR$Trait[1], '_relations.pdf'))
  psych::pairs.panels(subset(data_GR, select = c(Clim, Trait_mean,
                                                  Demog_rate_mean, Pop_mean, GR)),
                              ellipses = FALSE, hist.col = 'grey', lm = TRUE)
  mtext(paste0('Demographic rate is ', unique(data_GR$Demog_rate)), side = 3,
        line = 3)
  mtext(paste0('Trait is ', unique(data_GR$Trait_Categ_det)), side = 1,
        line = 4)
  dev.off()

  ## now call a function fitting a model (depending on the options:
  ## - autocor / no, DD/no, weights /no)

  SEM_mod <- fit_mod(biol_data = data_GR, ID = ID,
                     DD = DD,
                     weights = weights,
                     correlation = correlation) ## there is another problem with gls
  ## it does not return the object if called from within psem()


  ## then get output from the results of fitting the function
  res_SEM <- get_res_SEM(SEM_mod)

  res_SEM$DD <- DD
  res_SEM$weights <- weights
  res_SEM$corr <- correlation

  # write output
  saveRDS(object = res_SEM,
          file = paste0('./', out_SEM, '/', res_SEM$ID[1], '_',
                        res_SEM$Species[1], '_', res_SEM$Location[1],
                        '_', res_SEM$Trait[1], '_',
                        res_SEM$Demog_rate[1], '_ResultsSEM',  '.RDS'))
  return(res_SEM)
}
