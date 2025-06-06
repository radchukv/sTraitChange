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
#' @param all_Relations Vector specifying the names of the relations from the SEMs that will
#' be used as response variables in meta-analyses to be fitted.
#' @inheritParams fit_mod
#' @inheritParams fit_meta_phylo
#'
#' @export
#' @importFrom magrittr "%>%"
#'
#' @return Returns a tibble that includes one data frame and one tibble. A data frame
#' contains the estimated effect sizes, their standard errors, lower and upper confidence
#' intervals, their significance and the AIC for both the mixed-effects model fitted while
#' accounting and not for phylogenetic relatedness. This dataframe also has the column
#' 'Relation' specifying for which relation the data were analyzed
#' (e.g. Demog_rate_mean<-Pop_mean', for more details see \code{\link{fit_meta_phylo}}).
#' A tibble contains three lists and one character variable. The first list contains estimated global
#' effect sizes, their standard errors, their significance and the AIC per each fitted mixed-effects model
#' defined by the 'Relation', as returned by \code{\link{fit_meta_phylo}}. The second list
#' contains data frames with the effect sizes and standard errors for each study.
#' The third list contains the heterogeneity metrics extracted from each meta-analysis.
#' Each data frame corresponds to the subset of the data per each 'Relation' type.
#' A character 'names' specifies the 'Relation' type.
#'
#' @examples
#' # fit the models for: Trait_mean<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'
#' # rename the spp to have same names as in the phylogeny (nomenclature differs for some spp)
#' dataPaths <- dataPaths %>%
#'   dplyr::mutate(Sp_phylo = dplyr::case_when(
#'   Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
#'   Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
#'   Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
#'   Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
#'   Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
#'   Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
#'   TRUE ~ Species))
#'
#' # add underscore - same format as in phylogeny
#' dataPaths$Species <- unlist(lapply(1:nrow(dataPaths), FUN = function(x){
#' binary <- strsplit(as.character(dataPaths$Sp_phylo[x]), " ")
#' Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")
#' }))
#'
#' dataPaths$Sp_phylo <- dataPaths$Species
#' meta_Phen_Cov <- fit_all_meta(data_MA = dataPaths,
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
#' meta_Phen_Cov
#'
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
                         A,
                         des.matrix = "treatm.contrasts",
                         all_Relations = c('Demog_rate_mean<-det_Clim', 'Demog_rate_mean<-Pop_mean',
                                           'Demog_rate_mean<-Trait_mean', 'GR<-Demog_rate_mean',
                                           'GR<-det_Clim', 'GR<-Pop_mean',
                                           'Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim',
                                           'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim',
                                           'Trait_mean<-det_Clim', 'Ind_GR<-Pop_mean',
                                           'Tot_GR<-Pop_mean')){

  if(! is.null(Demog_rate)){
  subs_data <- droplevels(data_MA[data_MA$Demog_rate_Categ == Demog_rate & data_MA$Trait_Categ == Trait_categ,])
  } else {
    subs_data <- droplevels(data_MA[data_MA$Trait_Categ == Trait_categ,])
  }


  stat_meta <- do.call('rbind', lapply(1:length(all_Relations), FUN = function(x){
    fit_meta_phylo(data_MA = subs_data, Type_EfS = all_Relations[x], Cov_fact = Cov_fact,
             COV = COV, DD = DD, A = A, simpleSEM = simpleSEM,
             Trait = Trait, des.matrix = des.matrix)[, c('data','data_EfS', 'heter_mod')]}))


  rel_realized <- lapply(1:length(stat_meta$data_EfS), FUN = function(x){
    unique(stat_meta$data_EfS[[x]]$Relation)
  })


  stat_meta$names <- unlist(rel_realized)

  names(stat_meta$data) <- unlist(rel_realized)

  if (requireNamespace("data.table", quietly = TRUE)) {
  meta_res <- data.table::rbindlist(stat_meta$data)
  } else {
    message("to run fit_all_meta(), you first must run install.packages('data.table')!")
  }
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
    return(tibble::tibble(meta_res = list(meta_res),
                          data_meta = list(stat_meta)))

}
