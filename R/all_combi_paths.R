#' Use bootstrap to obtain the coefficients and SEs for all combined pathways
#'
#' \code{all_combi_paths} returns the medians and their SE estimates for
#' all combined pathways present in SEM, using the bootstrap procedure
#'
#' @param data_for_MA Data frame containing a subset of data_MA (see
#' \code{fit_meta_phylo}) in a wide format, i.e. single column per pathway.
#' This data frame contains for each pathway for each study the
#' effect size estimates along with their SEs obtained with the SEM
#' analyses. Besides that, some meta-data on each study is included
#' (e.g. study species, study location, study ID).
#' @inheritParams fit_mod
#'
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"
#'
#' @return Returns a data frame, which is an extended version of the input data frame,
#' that additionally contains estimated medians, SEs, as well as lower and upper CI
#' for the combined pathways.
#'
#' @examples
#' # prepare the data to fit the model to extract indirect path
#' dataPaths_sp <- dataPaths %>%
#'                   dplyr::mutate(Species = dplyr::case_when(
#'                          Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
#'                          Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
#'                          Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
#'                          Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
#'                          Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
#'                          Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
#'                          TRUE ~ Species))
#'
#' dataPaths_sp$Species <- unlist(lapply(1:nrow(dataPaths_sp), FUN = function(x){
#'   binary <- strsplit(as.character(dataPaths_sp$Species[x]), " ")
#'   Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")}))
#' dataPaths_sp$Sp_phylo <- dataPaths_sp$Species
#'
#' # fit the models for: Trait_mean<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'
#' meta_Phen_Cov <- fit_all_meta(data_MA = dataPaths_sp,
#'                               Demog_rate = NULL,
#'                               Trait_categ = 'Phenological',
#'                               Clim = 'Temperature',
#'                               Cov_fact = 'WeathQ',
#'                               COV = 'Pvalue',
#'                               sel = 'Temp_Phen_Cov',
#'                               folder_name = NULL,
#'                               colr = c('black', 'red'),
#'                               DD = 'n_effectGR',
#'                               simpleSEM = TRUE,
#'                               A = phyloMat,
#'                               all_Relations = c('Trait_mean<-det_Clim',
#'                               'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'))
#' two <- dplyr::bind_rows(lapply(list('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'),
#' function(x){extr_coefs(obj = meta_Phen_Cov, Type_EfS = x)}))
#' # fixing inconsistent variable naming
#' Coefs_subst <- dataPaths_sp %>%
#' dplyr::rename(.,  SError = Std.Error)
#' two$P.Value = rep(NA, nrow(two))
#'
#' Coefs_subst <- Coefs_subst %>%
#'     dplyr::select(., -c(names(Coefs_subst)[! names(Coefs_subst) %in% names(two)]))
#' two <- two %>%
#'     dplyr::select(., -c(names(two)[! names(two) %in% names(Coefs_subst)]))
#' Coef_all <- rbind(Coefs_subst, two) %>%
#' dplyr::filter(Trait_Categ == 'Phenological')
#'  forTrans <- subset(Coef_all, select = c(Estimate,  SError, Relation,
#'  Species, Location, ID))
#'  met_wide <- forTrans %>%
#'      tidyr::gather(variable, value, -(Relation:ID)) %>%
#'      tidyr::unite(temp, Relation, variable, sep = '/') %>%
#'      tidyr::spread(temp, value)
#'  test_simSEM <- all_combi_paths(data_for_MA = met_wide, DD = 'n_effectGR',
#'                                 simpleSEM = TRUE, Trait = FALSE)
#'
all_combi_paths <- function(data_for_MA, DD = 'n_effectDGR', simpleSEM = FALSE,
                            Trait = FALSE){

  if(simpleSEM){
    Ind_GR.df <- purrr::pmap_dfr(list(x = data_for_MA$`Trait_mean<-det_Clim/Estimate`,
                                      y = data_for_MA$`GR<-Trait_mean/Estimate`,
                                      x.se = data_for_MA$`Trait_mean<-det_Clim/SError`,
                                      y.se = data_for_MA$`GR<-Trait_mean/SError`),
                                 ind_path) %>%
      dplyr::rename(`Ind_GR<-det_Clim/Estimate` = .data$Median,
                    `Ind_GR<-det_Clim/SError` = .data$SE,
                    `Ind_GR<-det_Clim/lCI` = .data$lCI,
                    `Ind_GR<-det_Clim/uCI` = .data$uCI)

    Tot_GR.df <- purrr::pmap_dfr(list(direct = data_for_MA$`GR<-det_Clim/Estimate`,
                                      indir = Ind_GR.df$`Ind_GR<-det_Clim/Estimate`,
                                      direct.se = data_for_MA$`GR<-det_Clim/SError`,
                                      indir.se = Ind_GR.df$`Ind_GR<-det_Clim/SError`),
                                 tot_path)  %>%
      dplyr::rename(`Tot_GR<-det_Clim/Estimate` = .data$Median,
                    `Tot_GR<-det_Clim/SError` = .data$SE,
                    `Tot_GR<-det_Clim/lCI` = .data$lCI,
                    `Tot_GR<-det_Clim/uCI` = .data$uCI)
    met_wide <- cbind(data_for_MA, Ind_GR.df, Tot_GR.df)

  }else {
    if (! Trait){  ## Indirect path from climate to GR is different if the effect of trait on GR is also included
      ## in that case it is a sum of 2 paths: clim_Tr*Tr_DR*DR_GR and Tr_GR
    Ind_GR.df <- purrr::pmap_dfr(list(x = data_for_MA$`Trait_mean<-det_Clim/Estimate`,
                                      y = data_for_MA$`Demog_rate_mean<-Trait_mean/Estimate`,
                                      z = data_for_MA$`GR<-Demog_rate_mean/Estimate`,
                                      x.se = data_for_MA$`Trait_mean<-det_Clim/SError`,
                                      y.se = data_for_MA$`Demog_rate_mean<-Trait_mean/SError`,
                                      z.se = data_for_MA$`GR<-Demog_rate_mean/SError`),
                                 ind_path) %>%
      dplyr::rename(`Ind_GR<-det_Clim/Estimate` = .data$Median,
                    `Ind_GR<-det_Clim/SError` = .data$SE,
                    `Ind_GR<-det_Clim/lCI` = .data$lCI,
                    `Ind_GR<-det_Clim/uCI` = .data$uCI)
    } else {
      Ind_GR.df <- purrr::pmap_dfr(list(x = data_for_MA$`Trait_mean<-det_Clim/Estimate`,
                                        y = data_for_MA$`Demog_rate_mean<-Trait_mean/Estimate`,
                                        z = data_for_MA$`GR<-Demog_rate_mean/Estimate`,
                                        omega = data_for_MA$`GR<-Trait_mean/Estimate`,
                                        x.se = data_for_MA$`Trait_mean<-det_Clim/SError`,
                                        y.se = data_for_MA$`Demog_rate_mean<-Trait_mean/SError`,
                                        z.se = data_for_MA$`GR<-Demog_rate_mean/SError`,
                                        omega.se = data_for_MA$`GR<-Trait_mean/SError`),
                                   ind_path) %>%
        dplyr::rename(`Ind_GR<-det_Clim/Estimate` = .data$Median,
                      `Ind_GR<-det_Clim/SError` = .data$SE,
                      `Ind_GR<-det_Clim/lCI` = .data$lCI,
                      `Ind_GR<-det_Clim/uCI` = .data$uCI)

    }

    Ind_DemRate.df <- purrr::pmap_dfr(list(x = data_for_MA$`Trait_mean<-det_Clim/Estimate`,
                                           y = data_for_MA$`Demog_rate_mean<-Trait_mean/Estimate`,
                                           x.se = data_for_MA$`Trait_mean<-det_Clim/SError`,
                                           y.se = data_for_MA$`Demog_rate_mean<-Trait_mean/SError`),
                                      ind_path) %>%
      dplyr::rename(`Ind_DemRate<-det_Clim/Estimate` = .data$Median,
                    `Ind_DemRate<-det_Clim/SError` = .data$SE,
                    `Ind_DemRate<-det_Clim/lCI` = .data$lCI,
                    `Ind_DemRate<-det_Clim/uCI` = .data$uCI)

    Tot_DemRate.df <- purrr::pmap_dfr(list(direct = data_for_MA$`Demog_rate_mean<-det_Clim/Estimate`,
                                           indir = Ind_DemRate.df$`Ind_DemRate<-det_Clim/Estimate`,
                                           direct.se = data_for_MA$`Demog_rate_mean<-det_Clim/SError`,
                                           indir.se = Ind_DemRate.df$`Ind_DemRate<-det_Clim/SError`),
                                      tot_path) %>%
      dplyr::rename(`Tot_DemRate<-det_Clim/Estimate` = .data$Median,
                    `Tot_DemRate<-det_Clim/SError` = .data$SE,
                    `Tot_DemRate<-det_Clim/lCI` = .data$lCI,
                    `Tot_DemRate<-det_Clim/uCI` = .data$uCI)

    Tot_GR.df <- purrr::pmap_dfr(list(direct = data_for_MA$`GR<-det_Clim/Estimate`,
                                      indir = Ind_GR.df$`Ind_GR<-det_Clim/Estimate`,
                                      ClDem = data_for_MA$`Demog_rate_mean<-det_Clim/Estimate`,
                                      DemGR = data_for_MA$`GR<-Demog_rate_mean/Estimate`,
                                      direct.se = data_for_MA$`GR<-det_Clim/SError`,
                                      indir.se = Ind_GR.df$`Ind_GR<-det_Clim/SError`,
                                      ClDem.se = data_for_MA$`Demog_rate_mean<-det_Clim/SError`,
                                      DemGR.se = data_for_MA$`GR<-Demog_rate_mean/SError`),
                                 tot_path)  %>%
      dplyr::rename(`Tot_GR<-det_Clim/Estimate` = .data$Median,
                    `Tot_GR<-det_Clim/SError` = .data$SE,
                    `Tot_GR<-det_Clim/lCI` = .data$lCI,
                    `Tot_GR<-det_Clim/uCI` = .data$uCI)
    ## for DD acting on both DR and GR
    if (DD == 'n_effectDGR'){
      Ind_DD.df <- purrr::pmap_dfr(list(x = data_for_MA$`Demog_rate_mean<-Pop_mean/Estimate`,
                                        y = data_for_MA$`GR<-Demog_rate_mean/Estimate`,
                                        x.se = data_for_MA$`Demog_rate_mean<-Pop_mean/SError`,
                                        y.se = data_for_MA$`GR<-Demog_rate_mean/SError`),
                                   ind_path) %>%
        dplyr::rename(`Ind_GR<-Pop_mean/Estimate` = .data$Median,
                      `Ind_GR<-Pop_mean/SError` = .data$SE,
                      `Ind_GR<-Pop_mean/lCI` = .data$lCI,
                      `Ind_GR<-Pop_mean/uCI` = .data$uCI)

      Tot_DD.df <- purrr::pmap_dfr(list(direct = data_for_MA$`GR<-Pop_mean/Estimate`,
                                        indir = Ind_DD.df$`Ind_GR<-Pop_mean/Estimate`,
                                        direct.se = data_for_MA$`GR<-Pop_mean/SError`,
                                        indir.se = Ind_DD.df$`Ind_GR<-Pop_mean/SError`),
                                   tot_path) %>%
        dplyr::rename(`Tot_GR<-Pop_mean/Estimate` = .data$Median,
                      `Tot_GR<-Pop_mean/SError` = .data$SE,
                      `Tot_GR<-Pop_mean/lCI` = .data$lCI,
                      `Tot_GR<-Pop_mean/uCI` = .data$uCI)
      met_wide <- cbind(data_for_MA, Ind_GR.df, Ind_DemRate.df, Tot_GR.df, Tot_DemRate.df, Ind_DD.df, Tot_DD.df)

    }
    else{
      met_wide <- cbind(data_for_MA, Ind_GR.df, Ind_DemRate.df, Tot_GR.df, Tot_DemRate.df)
    }
  }

  return(met_wide)
}
