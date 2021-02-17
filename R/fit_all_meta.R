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
#' to subset the data. If NULL no subsetting by demographic rate is performed.
#' @param Trait_categ Character specifying the level of the trait on which to subset the data.
#' @param Clim Character specifying the level of the climatic variable on which to subset
#' the data.
#' @param sel Character specifying the name to be included in the .pdf name, which indicates
#' the levels of demographic rate and trait for which the subsetting was done.
#' @param folder_name Character specifyng the path to the directory in which the results will be saved.
#' @param colr Vector specifying the colours to be used for the data points. The length of the
#' vector should correspond to the number of the levels in the categorical explanatory variable
#' included in the meta-analytical model, or should be one (if the single global effect size
#' across all studies is to be plotted).
#' @inheritParams fit_mod
#' @inheritParams fit_meta
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
#'                                Demog_rate = 'survival',
#'                                Trait_categ = 'phenological',
#'                                Clim = 'temperature',
#'                                COV = NULL,
#'                                Covar = NULL,
#'                                sel = 'Phen_Surv',
#'                                folder_name = './output_overall/',
#'                                colr = c('black'))
#' meta_Phen_Surv
#' meta_Phen_Surv_byCont <- fit_all_meta(data_MA = Coefs_Aut,
#'                                       Demog_rate = 'survival',
#'                                       Trait_categ = 'phenological',
#'                                       Clim = 'temperature',
#'                                       COV = 'Pvalue',
#'                                       Covar = 'Continent',
#'                                       sel = 'Phen_Surv',
#'                                       folder_name = './output_overall/',
#'                                       colr = c('black', 'green', 'blue', 'red'))
fit_all_meta <- function(data_MA,
                         Demog_rate = 'Survival',
                         Trait_categ = 'Phenological',
                         Clim = 'Temperature',
                         Cov_fact = NULL,
                         COV = NULL,
                         sel = 'Phen_Surv',
                         folder_name = NULL,
                         colr = c('black'),
                         DD = 'n_effectDGR',
                         simpleSEM = FALSE,
                         Trait = FALSE,
                         all_Relations = c('Demog_rate_mean<-det_Clim', 'Demog_rate_mean<-Pop_mean',
                                           'Demog_rate_mean<-Trait_mean', 'GR<-Demog_rate_mean',
                                           'GR<-det_Clim', 'GR<-Pop_mean',
                                           'Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim',
                                           'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim',
                                           'Trait_mean<-det_Clim', 'Ind_GR<-Pop_mean',
                                           'Tot_GR<-Pop_mean')){

  ### subs_data <- droplevels(base::subset(data_MA, Demog_rate_Categ == Demog_rate & Trait_Categ == Trait_categ))  ## weird, why this is not working as intended???

  if(! is.null(Demog_rate)){
  subs_data <- droplevels(data_MA[data_MA$Demog_rate_Categ == Demog_rate & data_MA$Trait_Categ == Trait_categ,])
  } else {
    subs_data <- droplevels(data_MA[data_MA$Trait_Categ == Trait_categ,])
  }


  stat_meta <- do.call('rbind', lapply(1:length(all_Relations), FUN = function(x){
    fit_meta(data_MA = subs_data, Type_EfS = all_Relations[x], Cov_fact = Cov_fact,
             COV = COV, DD = DD, simpleSEM = simpleSEM, Trait = Trait)[, c('data','data_EfS', 'data_R2')]}))


  rel_realized <- lapply(1:length(stat_meta$data_EfS), FUN = function(x){
    unique(stat_meta$data_EfS[[x]]$Relation)
  })


  stat_meta$names <- unlist(rel_realized)

  names(stat_meta$data) <- unlist(rel_realized)
  meta_res <- data.table::rbindlist(stat_meta$data)
  if(is.null(COV)){
    if(! is.null(Cov_fact)){
      num_lev <- length(unique(meta_res$Levels_Covar))
      meta_res$Relation <- rep(unlist(rel_realized), each = num_lev)
    } else {
      meta_res$Relation <- unlist(rel_realized)
    }
  }  else{
    if(! is.null(Cov_fact)){
        num_lev <- length(unique(meta_res$Levels_Covar))

      meta_res$Relation <- rep(unlist(rel_realized), each = num_lev)
    } else {
      if(length(grep('\\+', COV)) == 1){
        meta_res$Relation <- rep(unlist(rel_realized), each = 3)
      } else {
        meta_res$Relation <- rep(unlist(rel_realized), each = 2)
      }
    }

  }

  ## a plot per each model
  lapply(1:length(unlist(rel_realized)), FUN = function(x){
    xlab <- plot_lab_name(Relation = stat_meta$names[[x]],
                          Cov_fact = Cov_fact,
                          Demog_rate = Demog_rate,
                          Trait_categ = Trait_categ,
                          Clim = Clim)$xlab
    coef <- plot_lab_name(Relation = stat_meta$names[[x]],
                          Cov_fact = Cov_fact,
                          Demog_rate = Demog_rate,
                          Trait_categ = Trait_categ,
                          Clim = Clim)$pref_pdf
    if(! is.null(folder_name)){
    if(is.null(Cov_fact)){
    plot_forest(data_ES = stat_meta$data_EfS[[x]],
                data_globES = stat_meta$data[[x]],
                Cov_fact = Cov_fact, xlab = xlab,
                colr = c('black'), COV = COV,
                pdf_basename = paste0(folder_name, sel, '_Coef_', coef),
                mar = c(4, 10, 2, 2),
                labels_ES = TRUE)
    }else {
      plot_forest(data_ES = stat_meta$data_EfS[[x]],
                  data_globES = stat_meta$data[[x]],
                  Cov_fact = Cov_fact, xlab = xlab,
                  colr = colr, COV = COV,
                  pdf_basename = paste0(folder_name, sel, '_Coef_', coef),
                  mar = c(4, 10, 2, 2),
                  labels_ES = TRUE)
    }
    }

  })

  ## getting prop. contributions
  if(any(all_Relations %in% c('Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim',
                          'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim',
                          'Ind_GR<-Pop_mean', 'Tot_GR<-Pop_mean'))){
    prop_data <- fit_meta(data_MA = subs_data, Type_EfS =
               all_Relations[all_Relations %in% c('Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim',
                                                  'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim',
                                                  'Ind_GR<-Pop_mean', 'Tot_GR<-Pop_mean')][1],
             Cov_fact = Cov_fact, COV = COV, DD = DD,
             simpleSEM = simpleSEM, Trait = Trait)[, 'prop_data']

    return(tibble::tibble(meta_res = list(meta_res),
                          data_meta = list(stat_meta),
                          prop_data = list(prop_data[[1]][[1]])))
  }
  else{
    return(tibble::tibble(meta_res = list(meta_res),
                          data_meta = list(stat_meta)))
  }


}
