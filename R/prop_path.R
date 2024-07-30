#' Obtain proportional contributions of each pathway to the total pathway
#'
#' \code{prop_path} estimates proportional contributions of each path (direct,
#' indirect and others if present) to the two total pathways: one from climate
#' to demographic rate and the other one from climate to population growth rate.
#' This function is called from within a function \code{\link{fit_meta}}.
#'
#' @param data_MA Data frame containing, for each study, the effect size
#' estimates for each pathway from the SEM analyses. This data frame also
#' contains meda-data (e.g. sutdy species, study location, continent, life
#' history traits of the species), that are needed to fit the mixed-effect
#' model. See also \code{\link{fit_meta}}
#' @param data Data frame containing, for each study, the effect size estimates
#' of all pathways, also combined ones. This dataframe is returned as an
#' intermediate step of the function \code{\link{fit_meta}}.
#' @inheritParams fit_mod
#'
#' @return Returns a data frame that in addition to the meta-data for each study
#' also contains proportional contributions of each pathway (direct, indirect,
#' and others, if present) to the total pathway from climate to demographic
#' rate (if \code{simpleSEM = FALSE}) and from climate to the population growth
#' rate.
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_forSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#' check_noCovar <- fit_meta(data_MA = Coefs_Aut, Type_EfS = 'Trait_mean<-det_Clim',
#'                           Covar = NULL)
#' check_TraitCateg <- fit_meta(data_MA = Coefs_Aut, Type_EfS = 'Trait_mean<-det_Clim',
#'                              Covar = 'Trait_Categ')
#' check_TraitCateg
prop_path <- function(data, data_MA, DD = 'n_effectDGR',
                      simpleSEM = FALSE){
  if(DD == 'n_effectDGR'){
  data_prop <- data %>%
    dplyr::mutate(abs_tot_Dem = abs(`Ind_DemRate<-det_Clim/Estimate`) +
                    abs(`Demog_rate_mean<-det_Clim/Estimate`),
                  abs_tot_GR = abs(`Ind_GR<-det_Clim/Estimate`) + abs(`GR<-det_Clim/Estimate`) +
                    abs(`Demog_rate_mean<-det_Clim/Estimate` * `GR<-Demog_rate_mean/Estimate`),
                  abs_tot_DD = abs(`Ind_GR<-Pop_mean/Estimate`) + abs(`GR<-Pop_mean/Estimate`),
                  prop_ind_Dem = abs(`Ind_DemRate<-det_Clim/Estimate`) / abs_tot_Dem,
                  prop_dir_Dem =  abs(`Demog_rate_mean<-det_Clim/Estimate`) / abs_tot_Dem,
                  prop_ind_GR = abs(`Ind_GR<-det_Clim/Estimate`) / abs_tot_GR,
                  prop_dir_GR = abs(`GR<-det_Clim/Estimate`) / abs_tot_GR,
                  prop_DemMed_GR = abs(`Demog_rate_mean<-det_Clim/Estimate` *
                                         `GR<-Demog_rate_mean/Estimate`) / abs_tot_GR,
                  prop_ind_DD = abs(`Ind_GR<-Pop_mean/Estimate`) / abs_tot_DD,
                  prop_dir_DD = abs(`GR<-Pop_mean/Estimate`) / abs_tot_DD)
  } else{
    if (simpleSEM){
      data_prop <- data %>%
        dplyr::mutate(abs_tot_GR = abs(`Ind_GR<-det_Clim/Estimate`) + abs(`GR<-det_Clim/Estimate`),
                      prop_ind_GR = abs(`Ind_GR<-det_Clim/Estimate`) / abs_tot_GR,
                      prop_dir_GR = abs(`GR<-det_Clim/Estimate`) / abs_tot_GR)
    }
        else {
          data_prop <- data %>%
            dplyr::mutate(abs_tot_Dem = abs(`Ind_DemRate<-det_Clim/Estimate`) +
                            abs(`Demog_rate_mean<-det_Clim/Estimate`),
                          abs_tot_GR = abs(`Ind_GR<-det_Clim/Estimate`) + abs(`GR<-det_Clim/Estimate`) +
                            abs(`Demog_rate_mean<-det_Clim/Estimate` * `GR<-Demog_rate_mean/Estimate`),
                          prop_ind_Dem = abs(`Ind_DemRate<-det_Clim/Estimate`) / abs_tot_Dem,
                          prop_dir_Dem =  abs(`Demog_rate_mean<-det_Clim/Estimate`) / abs_tot_Dem,
                          prop_ind_GR = abs(`Ind_GR<-det_Clim/Estimate`) / abs_tot_GR,
                          prop_dir_GR = abs(`GR<-det_Clim/Estimate`) / abs_tot_GR,
                          prop_DemMed_GR = abs(`Demog_rate_mean<-det_Clim/Estimate` *
                                                 `GR<-Demog_rate_mean/Estimate`) / abs_tot_GR)
        }

  }




  subs_merge <- droplevels(data_MA %>%
                             dplyr::distinct(., ID, Country, Continent,
                                             Longitude, Latitude, Taxon,
                                             BirdType,
                                             Trait_Categ, Trait,
                                             Demog_rate_Categ, Demog_rate,
                                             Count, Nyears, WinDur, deltaAIC,
                                             .keep_all = T) %>%
                             subset(.,
                                    select = c(ID, Country, Continent,
                                               Longitude, Latitude, Taxon,
                                               BirdType, Trait_Categ, Trait,
                                               Demog_rate_Categ, Demog_rate,
                                               Count, Nyears, WinDur, deltaAIC,
                                               Pvalue)))
  data_prop <-  merge(data_prop, subs_merge, by = c('ID'))

  return(data_prop)
}
