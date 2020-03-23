#' Fit all meta-analytical models for each relation in the SEM,
#' for a specified subselection of demographic rate and trait
#'
#' \code{fit_all_meta} fits all mixed-effects meta-analytical models (per each relation
#' in the SEM), to extract global effect sizes across the studies in a specified subset
#' of data as defined by a demographic rate and a trait
#'
#' @param data_MA Data frame containing, for each study, the effect size estimates
#' for each pathway from the SEM analyses. This data frame also contains meda-data
#' (e.g. sutdy species, study location, continent, life history traits of the species),
#' that are needed to fit the mixed-effect model.
#' @param Demog_rate Character specifying the level of the demographic rate on which
#' to subset the data.
#' @param Trait_categ Character specifying the level of the trait on which to subset the data.
#' @param Clim Character specifying the level of the climatic variable on which to subset
#' the data.
#' @param tab Categorical specifying the name of the categorical variable for which
#'  to produce the frequency table.
#' @param Covar Categorical specifying the name of the categorical variable to be included
#' as fixed-effect covariate in the meta-analysis. Defaults to NULL, in which case no
#' categorical variables are included and the overall global effect size is estimated.
#' @param sel Character specifying the name to be included in the .pdf name, which indicates
#' the levels of demographic rate and trait for which the subsetting was done.
#' @param folder_name Character specifyng the path to the directory in which the results will be saved.
#'
#' @export
#'
#' @return Returns a tibble that includes one data frame, one tibble and a table. A data frame
#' contains the estimated effect sizes, their standard errors, their significance and the AIC per each
#' fitted mixed-effects model, including also the column 'Relation' specifying for which relation
#' the data were analyzed (e.g. Demog_rate_mean<-Pop_mean', for more details see \code{\link{fit_meta}}).
#' A tibble contains two lists and one character variable. The first list contains estimated global
#' effect sizes, their standard errors, their significance and the AIC per each fitted mixed-effects model
#' defined by the 'Relation', as returned by \code{\link{fit_meta}}. The second list contains data frames
#' with the effect sizes and standard errors for each study. Each data frame corresponds to the subset
#' of the data per each 'Relation' type. A character 'names' specifies the 'Relation' type. A table
#' contains frequencies for the levels of the specified categorical variable, supplied via parameter 'tab'.
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_forSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#' meta_Phen_Surv <- fit_all_meta(data_MA = Coefs_Aut,
#'                                Demog_rate = 'Survival',
#'                                Trait_categ = 'Phenological',
#'                                tab = 'Taxon',
#'                                Covar = NULL)
#' meta_Phen_Surv

fit_all_meta <- function(data_MA,
                         Demog_rate = 'survival',
                         Trait_categ = 'phenological',
                         Clim = 'temperature',
                         tab = 'Taxon',
                         Covar = NULL,
                         sel = 'Phen_Surv',
                         folder_name = './output_overall/'){

  ### subs_data <- droplevels(base::subset(data_MA, Demog_rate_Categ1 == Demog_rate & Trait_Categ == Trait_categ))  ## weird, why this is not working as intended???

  subs_data <- droplevels(data_MA[data_MA$Demog_rate_Categ1 == Demog_rate & data_MA$Trait_Categ == Trait_categ,])
  all_Relations <- c('Demog_rate_mean<-det_Clim', 'Demog_rate_mean<-Pop_mean', 'Demog_rate_mean<-Trait_mean',
                     'GR<-Demog_rate_mean',  'GR<-det_Clim', 'GR<-Pop_mean', 'Ind_DemRate<-det_Clim',
                     'Ind_GR<-det_Clim', 'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim', 'Trait_mean<-det_Clim')

  stat_meta <- do.call('rbind', lapply(all_Relations, FUN = function(x){
    fit_meta(data_MA = subs_data, Type_EfS = x, Covar = NULL)}))

  rel_realized <- lapply(1:length(stat_meta$data_EfS), FUN = function(x){
    unique(stat_meta$data_EfS[[x]]$Relation)
  })


  stat_meta$names <- unlist(rel_realized)
  names(stat_meta$data) <- unlist(rel_realized)
  meta_res <- data.table::rbindlist(stat_meta$data)
  meta_res$Relation <- unlist(rel_realized)
  data_fin <- stat_meta$data_EfS[[1]]  ## this may need ot be done cleaner or explain that I only use this subset for the descriptive stats on the factors

  tab_Taxtable <- table(data_fin[, tab])

  ## a plot per each model
  lapply(1:length(unlist(rel_realized)), FUN = function(x){
    xlab <- plot_lab_name(Relation = stat_meta$names[[x]],
                          Covar = Covar,
                          Demog_rate = Demog_rate,
                          Trait_categ = Trait_categ,
                          Clim = Clim)$xlab
    coef <- plot_lab_name(Relation = stat_meta$names[[x]],
                          Covar = Covar,
                          Demog_rate = Demog_rate,
                          Trait_categ = Trait_categ,
                          Clim = Clim)$pref_pdf
    plot_forest(data_ES = stat_meta$data_EfS[[x]],
                data_globES = stat_meta$data[[x]],
                Covar = NULL, xlab = xlab,
                colr = c('black'),
                pdf_basename = paste0(folder_name, sel, '_Coef_', coef),
                mar = c(4, 10, 2, 2),
                labels_ES = TRUE)

  })

  return(tibble::tibble(meta_res = list(meta_res),
                        data_meta = list(stat_meta),
                        prop_tab = list(tab_Taxtable)))
}
