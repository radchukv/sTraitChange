#' Fit a meta-analytical model to assess the global effects across the studies
#'
#' \code{fit_meta} fits mixed-effects meta-analytical model to extract a global
#' effect size across the studies
#'
#' @param data_MA Data frame containing, for each study, the effect size estimates
#' for each pathway from the SEM analyses. This data frame also contains meda-data
#' (e.g. sutdy species, study location, continent, life history traits of the species),
#' that are needed to fit the mixed-effect model.
#' @param Type_EfS Character specifying which effect size type to use for a meta-analysis.
#' These types of effect sizes reflect different pathways in the fitted SEM. Possible are:
#' 'Demog_rate_mean<-det_Clim', 'Demog_rate_mean<-Pop_mean', 'Demog_rate_mean<-Trait_mean',
#' 'GR<-Demog_rate_mean',  'GR<-det_Clim', 'GR<-Pop_mean', 'Ind_DemRate<-det_Clim',
#' 'Ind_GR<-det_Clim', 'Tot_DemRate<-det_Clim', 'Tot_GR<-Clim', and 'Trait_mean<-det_Clim'.
#' @param Covar Categorical specifying the name of the categorical variable to be included
#' as fixed-effect covariate in the meta-analysis. Defaults to NULL, in which case no
#' categorical variables are included and the overall global effect size is estimated.
#'
#' @export
#'
#' @return Returns a tibble. If no categorical explanatory variables are included in the model,
#' then this tibble includes two data frames: the data frame with the estimate of the global
#' effect size, its standard error, its significance and the AIC of the fitted model, and the
#' data frame with the effect sizes and their standard errors for each study.
#' If a categorical explanatory variable is included in the model, then the tibble
#' additionally contains the estimates of effect sizes for each level of the categorical variable,
#' as well as their standard errors, the p value indicating whether the categorical variable is
#' significant and the AIC of the model with categorical variable. The AIC and the p values are obtained
#' from the models fitted with ML, whereas the parameter estimates from the models fitted using REML.
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_forSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#' check_noCovar <- fit_meta(data_MA = Coefs_Aut, Type_EfS = 'Trait_mean<-det_Clim',
#'                           Covar = NULL)
#' check_TraitCateg <- fit_meta(data_MA = Coefs_Aut, Type_EfS = 'Trait_mean<-det_Clim',
#'                              Covar = 'Trait_Categ')
#' check_TraitCateg
fit_meta <- function(data_MA, Type_EfS = 'Trait_mean<-det_Clim',
                     Covar = NULL){
  ## calculating indirect effects, total effects and their SEs
  forTrans <- subset(data_MA, select = c(Estimate,  Std.Error, Relation, Species, Location, ID))
  forTrans <- forTrans %>%
    dplyr::rename(SError = Std.Error)

  met_wide <- forTrans %>%
    tidyr::gather(variable, value, -(Relation:ID)) %>%
    tidyr::unite(temp, Relation, variable, sep = '/') %>%
    tidyr::spread(temp, value)



  met_wide$`Ind_DemRate<-det_Clim/Estimate` <- met_wide$`Trait_mean<-det_Clim/Estimate` *
    met_wide$`Demog_rate_mean<-Trait_mean/Estimate`
  met_wide$`Ind_GR<-det_Clim/Estimate` <- met_wide$`Trait_mean<-det_Clim/Estimate` *
    met_wide$`Demog_rate_mean<-Trait_mean/Estimate` * met_wide$`GR<-Demog_rate_mean/Estimate`
  met_wide$`Tot_DemRate<-det_Clim/Estimate` <- met_wide$`Ind_DemRate<-det_Clim/Estimate` +
    met_wide$`Demog_rate_mean<-det_Clim/Estimate`
  met_wide$`Tot_GR<-Clim/Estimate` <- met_wide$`Ind_GR<-det_Clim/Estimate` + met_wide$`GR<-det_Clim/Estimate` +
    met_wide$`Demog_rate_mean<-det_Clim/Estimate` * met_wide$`GR<-Demog_rate_mean/Estimate`

  ## calculation of SE for those composite Ef sizes (like indirect effect and the total effect)
  ## using the formula from https://www.stata.com/statalist/archive/2005-12/msg00165.html
  ## implement later on as a funciton for a better use
  met_wide$`Ind_DemRate<-det_Clim/SError` <-  sqrt(met_wide$`Trait_mean<-det_Clim/Estimate`^2*met_wide$`Demog_rate_mean<-Trait_mean/SError`^2 +
                                                     met_wide$`Demog_rate_mean<-Trait_mean/Estimate`^2*met_wide$`Trait_mean<-det_Clim/SError`^2 +
                                                     met_wide$`Demog_rate_mean<-Trait_mean/SError`^2*met_wide$`Trait_mean<-det_Clim/SError`^2)  #a^2*V(b) + b^2*V(a) + V(a)*V(b)

  ## check this formula by calculating the SE also using the general formula- produces exactly the same result as above formula
  # met_wide$`Ind_Demrate<-Clim/SECheck` <- sqrt((met_wide$`Trait_mean<-det_Clim/SError`^2 + met_wide$`Trait_mean<-det_Clim/Estimate`^2)*
  #                                                (met_wide$`Demog_rate_mean<-Trait_mean/SError`^2 + met_wide$`Demog_rate_mean<-Trait_mean/Estimate`^2) -
  #                                                (met_wide$`Trait_mean<-det_Clim/Estimate`^2*met_wide$`Demog_rate_mean<-Trait_mean/Estimate`^2))

  ## I use this formula: https://stats.stackexchange.com/questions/52646/variance-of-product-of-multiple-random-variables
  met_wide$`Ind_GR<-det_Clim/SError` <- sqrt((met_wide$`Trait_mean<-det_Clim/SError`^2 + met_wide$`Trait_mean<-det_Clim/Estimate`^2)*
                                               (met_wide$`Demog_rate_mean<-Trait_mean/SError`^2 + met_wide$`Demog_rate_mean<-Trait_mean/Estimate`^2)*
                                               (met_wide$`GR<-Demog_rate_mean/SError`^2 + met_wide$`GR<-Demog_rate_mean/Estimate`^2) -
                                               (met_wide$`Trait_mean<-det_Clim/Estimate`^2 * met_wide$`Demog_rate_mean<-Trait_mean/Estimate`^2 *
                                                  met_wide$`GR<-Demog_rate_mean/Estimate`^2))

  ## the VAR of the sum of idnependent random vars is the sum of the variances
  met_wide$`Tot_DemRate<-det_Clim/SError` <- sqrt(met_wide$`Demog_rate_mean<-det_Clim/SError`^2 +
                                                    met_wide$`Ind_DemRate<-det_Clim/SError`^2)
  met_wide$`Tot_GR<-Clim/SError` <- sqrt(met_wide$`GR<-det_Clim/SError`^2 + met_wide$`Ind_GR<-det_Clim/SError`^2 +
                                           (met_wide$`Demog_rate_mean<-det_Clim/Estimate`^2*met_wide$`GR<-Demog_rate_mean/SError`^2 +
                                              met_wide$`GR<-Demog_rate_mean/Estimate`^2*met_wide$`Demog_rate_mean<-det_Clim/SError`^2 +
                                              met_wide$`GR<-Demog_rate_mean/SError`^2*met_wide$`Demog_rate_mean<-det_Clim/SError`^2))


  trans_allEfS <- met_wide %>%
    tidyr::gather(key, value, -c(Species:ID)) %>%
    tidyr::separate(., key, into = c('Relation', 'Metric'), sep = "/") %>%
    tidyr::spread(., Metric, value)

  subs_merge <- droplevels(subset(data_MA,
                                  select = c(ID, Country, Continent,
                                             Longitude, Latitude, Taxon,
                                             Trait_Categ, Trait, Demog_rate_Categ,
                                             Demog_rate_Categ1, Demog_rate, Count,
                                             Nyears, WinDur, deltaAIC, Pvalue)) %>%
                             dplyr::distinct(., ID, Country, Continent,
                                             Longitude, Latitude, Taxon,
                                             Trait_Categ, Trait, Demog_rate_Categ,
                                             Demog_rate_Categ1, Demog_rate, Count,
                                             Nyears, WinDur, deltaAIC, .keep_all = T))

  tot <- merge(trans_allEfS, subs_merge, by = 'ID')



  ## subset a specified effect size only
  subs_data <- subset(tot, Relation == Type_EfS)  ## here subset only by the type of effect size, the other subsetting should be done
  ## before fitting the data to this function(e.g. by taxon, trait categ etc.)


  ## if no categorical explanatories included
  mod_ML <- metafor::rma.mv(Estimate ~ 1, V = SError^2, #W = 1 / Pvalue,  ## works just as well as yi = Estimate
                            random = list(~ 1|Species, ~1|ID, ~1|Location), data = subs_data, method = 'ML')
  mod_REML <- metafor::rma.mv(yi = Estimate, V = SError^2, #W = 1 / Pvalue,
                              random = list(~ 1|Species, ~1|ID, ~1|Location), data = subs_data, method = 'REML')
  ## now trying to add another weight based on inverse of the p value - works, but has to be checked properly
  #mod_test_W <- metafor::rma.mv(yi = Estimate, V = SError^2, W = 1 / Pvalue,
  #                            random = list(~ 1|Species, ~1|ID, ~1|Location), data = subs_data, method = 'ML')
  ## this still has to be finished


  ## including the covariate
  if(! is.null(Covar)){
    formul <- paste0('Estimate ~ ', Covar, ' - 1')
    mod_ML_Cov <- metafor::rma.mv(stats::as.formula(formul), V = SError^2, #W = 1 / Pvalue,
                                  random = list(~ 1|Species, ~1|ID, ~1|Location), data = subs_data, method = 'ML')
    mod_REML_Cov <- metafor::rma.mv(stats::as.formula(formul), V = SError^2, #W = 1 / Pvalue,
                                    random = list(~ 1|Species, ~1|ID, ~1|Location), data = subs_data, method = 'ML')
    LRT_test <- metafor::anova.rma(mod_ML_Cov, mod_ML)
  }

  out_dat <- data.frame(EfS_across = as.numeric(mod_REML$beta), SE_across = mod_REML$se,
                        EfS_Low = mod_REML$ci.lb, EfS_Upper = mod_REML$ci.ub, pval_across =
                          mod_ML$pval, AIC_EfS_across = AIC(mod_ML))
  if(! is.null(Covar)){
    out_tib <- tibble::tibble(data = list(out_dat),
                              data_EfS = list(subs_data),
                              EfS_Covar = list(mod_REML_Cov$beta),
                              EfS_Covar_SE = list(mod_REML_Cov$se),
                              EfS_Covar_Low = list(mod_REML_Cov$ci.lb),
                              EfS_Covar_High = list(mod_REML_Cov$ci.ub),
                              pval_Covar = LRT_test$pval, AIC_EfS_Covar = AIC(mod_ML_Cov))
  }else{
    out_tib <- tibble::tibble(data = list(out_dat),
                              data_EfS = list(subs_data))
  }
  return(out_tib)
}
