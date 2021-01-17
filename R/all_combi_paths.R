#' Use bootstrap to obtain the coefficients and SEs for all combined pathways
#'
#' \code{all_combi_paths} returns the medians and their SE estimates for
#' all combined pathways present in SEM, using the bootstrap procedure
#'
#' @param data_for_MA Data frame containing a subset of data_MA (see
#' \code{fit_meta}) in a wide format, i.e. single column per pathway.
#' This data frame contains for each pathway for each study the
#' effect size estimates along with their SEs obtained with the SEM
#' analyses. Besides that, some meta-data on each study is included
#' (e.g. study species, study location, study ID).
#' @inheritParams fit_mod
#'
#' @export
#'
#' @return Returns a data frame, which is an extended version of the input data frame,
#' that addittionally contains estimates medians and SEs for the combined pathways.
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_fSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#' forTrans <- subset(Coefs_Aut, select = c(Estimate,  Std.Error, Relation, Species, Location, ID))
#' forTrans <- forTrans %>%
#'             dplyr::rename(SError = Std.Error)
#' met_wide <- forTrans %>%
#'             tidyr::gather(variable, value, -(Relation:ID)) %>%
#'             tidyr::unite(temp, Relation, variable, sep = '/') %>%
#'             tidyr::spread(temp, value)
#' test <- all_combi_paths(data = met_wide, DD = 'n_effectDGR', simpleSEM = FALSE)
#'
#' Coefs_Aut_simSEM <- readRDS(file = './output_fSEM_temp_simpleSEM/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#' forTrans_simSEM <- subset(Coefs_Aut_simSEM, select = c(Estimate,  Std.Error, Relation, Species, Location, ID))
#' forTrans_simSEM <- forTrans_simSEM %>%
#'                    dplyr::rename(SError = Std.Error)
#' met_wide_simSEM <- forTrans_simSEM %>%
#'                    tidyr::gather(variable, value, -(Relation:ID)) %>%
#'                    tidyr::unite(temp, Relation, variable, sep = '/') %>%
#'                    tidyr::spread(temp, value)
#' test_simSEM <- all_combi_paths(data = met_wide, DD = 'n_effectGR', simpleSEM = TRUE)
#'
all_combi_paths <- function(data_forMA, DD = 'n_effectDGR', simpleSEM = FALSE){

  if(simpleSEM){
    Ind_GR.df <- purrr::pmap_dfr(list(x = data_forMA$`Trait_mean<-det_Clim/Estimate`,
                                      y = data_forMA$`GR<-Trait_mean/Estimate`,
                                      x.se = data_forMA$`Trait_mean<-det_Clim/SError`,
                                      y.se = data_forMA$`Demog_rate_mean<-Trait_mean/SError`,
                                      z.se = data_forMA$`GR<-Trait_mean/SError`),
                                 ind_path) %>%
      dplyr::rename(., `Ind_GR<-det_Clim/Estimate` = Median,
                    `Ind_GR<-det_Clim/SError` = SE,
                    `Ind_GR<-det_Clim/lCI` = lCI,
                    `Ind_GR<-det_Clim/uCI` = uCI)

    Tot_GR.df <- purrr::pmap_dfr(list(direct = data_forMA$`GR<-det_Clim/Estimate`,
                                      indir = Ind_GR.df$`Ind_GR<-det_Clim/Estimate`,
                                      direct.se = data_forMA$`GR<-det_Clim/SError`,
                                      indir.se = Ind_GR.df$`Ind_GR<-det_Clim/SError`),
                                 tot_path)  %>%
      dplyr::rename(., `Tot_GR<-det_Clim/Estimate` = Median,
                    `Tot_GR<-det_Clim/SError` = SE,
                    `Tot_GR<-det_Clim/lCI` = lCI,
                    `Tot_GR<-det_Clim/uCI` = uCI)
    met_wide <- cbind(data_forMA, Ind_GR.df, Tot_GR.df)

  }else {
    Ind_GR.df <- purrr::pmap_dfr(list(x = data_forMA$`Trait_mean<-det_Clim/Estimate`,
                                      y = data_forMA$`Demog_rate_mean<-Trait_mean/Estimate`,
                                      z = data_forMA$`GR<-Demog_rate_mean/Estimate`,
                                      x.se = data_forMA$`Trait_mean<-det_Clim/SError`,
                                      y.se = data_forMA$`Demog_rate_mean<-Trait_mean/SError`,
                                      z.se = data_forMA$`GR<-Demog_rate_mean/SError`),
                                 ind_path) %>%
      dplyr::rename(., `Ind_GR<-det_Clim/Estimate` = Median,
                    `Ind_GR<-det_Clim/SError` = SE,
                    `Ind_GR<-det_Clim/lCI` = lCI,
                    `Ind_GR<-det_Clim/uCI` = uCI)

    Ind_DemRate.df <- purrr::pmap_dfr(list(x = data_forMA$`Trait_mean<-det_Clim/Estimate`,
                                           y = data_forMA$`Demog_rate_mean<-Trait_mean/Estimate`,
                                           x.se = data_forMA$`Trait_mean<-det_Clim/SError`,
                                           y.se = data_forMA$`Demog_rate_mean<-Trait_mean/SError`),
                                      ind_path) %>%
      dplyr::rename(., `Ind_DemRate<-det_Clim/Estimate` = Median,
                    `Ind_DemRate<-det_Clim/SError` = SE,
                    `Ind_DemRate<-det_Clim/lCI` = lCI,
                    `Ind_DemRate<-det_Clim/uCI` = uCI)

    Tot_DemRate.df <- purrr::pmap_dfr(list(direct = data_forMA$`Demog_rate_mean<-det_Clim/Estimate`,
                                           indir = Ind_DemRate.df$`Ind_DemRate<-det_Clim/Estimate`,
                                           direct.se = data_forMA$`Demog_rate_mean<-det_Clim/SError`,
                                           indir.se = Ind_DemRate.df$`Ind_DemRate<-det_Clim/SError`),
                                      tot_path) %>%
      dplyr::rename(., `Tot_DemRate<-det_Clim/Estimate` = Median,
                    `Tot_DemRate<-det_Clim/SError` = SE,
                    `Tot_DemRate<-det_Clim/lCI` = lCI,
                    `Tot_DemRate<-det_Clim/uCI` = uCI)

    Tot_GR.df <- purrr::pmap_dfr(list(direct = data_forMA$`GR<-det_Clim/Estimate`,
                                      indir = Ind_GR.df$`Ind_GR<-det_Clim/Estimate`,
                                      ClDem = data_forMA$`Demog_rate_mean<-det_Clim/Estimate`,
                                      DemGR = data_forMA$`GR<-Demog_rate_mean/Estimate`,
                                      direct.se = data_forMA$`GR<-det_Clim/SError`,
                                      indir.se = Ind_GR.df$`Ind_GR<-det_Clim/SError`,
                                      ClDem.se = data_forMA$`Demog_rate_mean<-det_Clim/SError`,
                                      DemGR.se = data_forMA$`GR<-Demog_rate_mean/SError`),
                                 tot_path)  %>%
      dplyr::rename(., `Tot_GR<-det_Clim/Estimate` = Median,
                    `Tot_GR<-det_Clim/SError` = SE,
                    `Tot_GR<-det_Clim/lCI` = lCI,
                    `Tot_GR<-det_Clim/uCI` = uCI)
    ## for DD acting on both DR and GR
    if (DD == 'n_effectDGR'){
      Ind_DD.df <- purrr::pmap_dfr(list(x = data_forMA$`Demog_rate_mean<-Pop_mean/Estimate`,
                                        y = data_forMA$`GR<-Demog_rate_mean/Estimate`,
                                        x.se = data_forMA$`Demog_rate_mean<-Pop_mean/SError`,
                                        y.se = data_forMA$`GR<-Demog_rate_mean/SError`),
                                   ind_path) %>%
        dplyr::rename(., `Ind_GR<-Pop_mean/Estimate` = Median,
                      `Ind_GR<-Pop_mean/SError` = SE,
                      `Ind_GR<-Pop_mean/lCI` = lCI,
                      `Ind_GR<-Pop_mean/uCI` = uCI)

      Tot_DD.df <- purrr::pmap_dfr(list(direct = data_forMA$`GR<-Pop_mean/Estimate`,
                                        indir = Ind_DD.df$`Ind_GR<-Pop_mean/Estimate`,
                                        direct.se = data_forMA$`GR<-Pop_mean/SError`,
                                        indir.se = Ind_DD.df$`Ind_GR<-Pop_mean/SError`),
                                   tot_path) %>%
        dplyr::rename(., `Tot_GR<-Pop_mean/Estimate` = Median,
                      `Tot_GR<-Pop_mean/SError` = SE,
                      `Tot_GR<-Pop_mean/lCI` = lCI,
                      `Tot_GR<-Pop_mean/uCI` = uCI)
      met_wide <- cbind(data_forMA, Ind_GR.df, Ind_DemRate.df, Tot_GR.df, Tot_DemRate.df, Ind_DD.df, Tot_DD.df)

    }
    else{
      met_wide <- cbind(data_forMA, Ind_GR.df, Ind_DemRate.df, Tot_GR.df, Tot_DemRate.df)
    }
  }

  return(met_wide)
}
