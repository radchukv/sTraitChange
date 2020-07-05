#' Extract results of the fitted structural equation model
#'
#' \code{get_res_SEM} extracts results of the fitted
#'  structural equation model
#'
#' @param mod_obj fitted psem object, returned with function
#' \code{\link{fit_mod}}
#'
#' @export
#'
#' @return A tibble with 12 columns that represent the metadata for a
#' given study ID and the results of SEM: Cstat, dTable, R2 of all
#' equations and path coefficients.
#'
#' @examples
#' biol_dat <- read.csv('./data-raw/Closed_Data_26_02.csv')
#' biol_stand <- convert_JulianDay(biol_data = biol_dat)
#' # example for ID = 1
#' subs <- droplevels(biol_stand[biol_stand$ID == 1, ])
#' test_fSEM <- readRDS(paste0('./output_forSEM_temp/', subs$ID[1], '_',
#'                     subs$Species[1], '_', subs$Location[1],
#'                     '_', subs$Trait[1], '_ForSEM',  '.RDS'))
#' dat <- test_fSEM$data_res[[1]]
#' ## impute the DR
#' dat <- impute_ma(data = dat, column = 'Demog_rate_mean')
#' dat <- impute_median(data = dat, column = 'Demog_rate_SE')
#' ## impute the Pop Size
#' dat <- impute_ma(data = dat, column = 'Pop_mean')
#' if(unique(subs$Count) == 'N'){
#' dat <- impute_median(data = dat, column = 'Pop_SE')
#' }
#' ## calculate GR
#' data_GR <- dat %>%
#' dplyr::mutate(., Pop_mean_lag = lag(Pop_mean)) %>%
#' dplyr::mutate(., GR = Pop_mean / Pop_mean_lag) %>%
#' dplyr::filter(., !is.na(GR))
#'
#' # fit the model
#' mod_SEM <- fit_mod(biol_data = data_GR, ID =1,
#'                 DD = FALSE, weights = FALSE,
#'                 correlation = FALSE)  ## does not work for corr = TRUE yet
#' res_SEM <- get_res_SEM(mod_SEM)

get_res_SEM <- function(mod_obj){

  summary_mod <- summary(mod_obj, .progressBar = F)
  dTable <- summary_mod$dTable
  Cstat <- summary_mod$Cstat
  R2_Relation <- piecewiseSEM::rsquared(mod_obj)
  coefs <- piecewiseSEM::coefs(mod_obj)

  AIC <- stats::AIC(mod_obj, aicc = TRUE)
  data <- mod_obj$data

  # add timeseries duration
  data <- data %>%
    dplyr::group_by(., ID) %>%
    dplyr::mutate(., NYears = dplyr::n()) %>%
    dplyr::ungroup()

  res <- tibble::tibble(ID = unique(data$ID),
                        Study_Authors = unique(data$Study_Authors),
                        Species = unique(data$Species),
                        Location = unique(data$Location),
                        Country = unique(data$Country),
                        Continent = unique(data$Continent),
                        Longitude = unique(data$Longitude),
                        Latitude = unique(data$Latitude),
                        Taxon = unique(data$Taxon),
                        BirdType = unique(data$BirdType),
                        Trait_Categ = unique(data$Trait_Categ),
                        Trait = unique(data$Trait),
                        Demog_rate_Categ = unique(data$Demog_rate_Categ),
                        Demog_rate = unique(data$Demog_rate),
                        Count = unique(data$Count),
                        Nyears = unique(data$NYears),
                        WinDur = unique(data$WinDur),
                        deltaAIC = unique(data$deltaAIC),
                        Pvalue = unique(data$Pvalue),
                        dTable = list(dTable),
                        Cstat = list(Cstat),
                        R2 = list(R2_Relation),
                        coefs = list(coefs),
                        AIC = AIC)

  return(res)
}
