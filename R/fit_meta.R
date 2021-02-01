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
#' 'Ind_GR<-det_Clim', 'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim', and 'Trait_mean<-det_Clim'.
#' @param Cov_fact A character specifying the name of the categorical variable to be included
#' as a fixed-effect covariate in the meta-analysis. Defaults to NULL, in which case no
#' categorical variables are included and the overall global effect size is estimated.
#' @param COV A character specifying the name (or names, separated by '+') of the continuous
#' variable to be included as a fixed-effect covariate in the meta-analysis. Defaults to NULl,
#' in which case no continuous covariates are included and the overall global effect size
#' is estimated.
#' @param optimize A character specifying the name of the optimizer to use. For more details
#' see argument 'control = ' for the available optimizers in \code{\link[metafor]{rma.mv}}.
#' Defaults to 'nlminb' because this optimizer was the most successful for model convergence
#' for our datasets.
#' @inheritParams fit_mod
#'
#' @export
#'
#' @return Returns a tibble. If no categorical explanatory variables are included in the model,
#' then this tibble includes four data frames: 1. 'data' - the data frame with the estimate
#' of the global effect size, its standard error, its significance, the AIC of the fitted
#' mixed-effects model and the variances explained by random effects;
#' 2. 'data_Efs' - the data frame with the effect sizes and their standard errors for each study;
#' 3. 'data_R2' - the data frame that includes in addition to the effect size and standard error
#' per each study also R2 of the fitted relations; and 4. 'prop_data' the data frame with
#' the proportional contribution of the direct and indirect paths to the total path. Additionally,
#' the column 'names' contains the names of the fitted relation for each mixed-effects model.
#' If a categorical explanatory variable is included in the model, then the tibble
#' additionally contains the dara frame with the estimates of effect sizes (and their SEs)
#' for each level of the categorical variable; and the columns with the p value
#' indicating whether the categorical variable is significant, the AIC of the mixed-effects model
#' with categorical variable and the variances explained by random effects.
#' The AIC and the p values are obtained from the models fitted with ML, whereas
#' the parameter estimates from the models fitted using REML.
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_forSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#' check_noCovar <- fit_meta(data_MA = Coefs_Aut, Type_EfS = 'Trait_mean<-det_Clim',
#'                           Cov_fact = NULL)
#' check_TraitCateg <- fit_meta(data_MA = Coefs_Aut, Type_EfS = 'Trait_mean<-det_Clim',
#'                              Cov_fact = 'Trait_Categ')
#' check_TraitCateg
fit_meta <- function(data_MA, Type_EfS = 'Trait_mean<-det_Clim',
                     Cov_fact = NULL, COV = NULL, optimize = 'uobyqa',
                     DD = 'n_effectDGR', simpleSEM = FALSE,
                     Trait = FALSE){
  ## calculating indirect effects, total effects and their SEs
  forTrans <- subset(data_MA, select = c(Estimate,  Std.Error, Relation, Species, Location, ID))
  forTrans <- forTrans %>%
    dplyr::rename(SError = Std.Error)

  met_wide <- forTrans %>%
    tidyr::gather(variable, value, -(Relation:ID)) %>%
    tidyr::unite(temp, Relation, variable, sep = '/') %>%
    tidyr::spread(temp, value)

  ## these should be called only if the Type_EfS is one of the
  ## indirect or total paths, to save time
  if(Type_EfS %in% c('Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim',
                     'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim',
                     'Ind_GR<-Pop_mean', 'Tot_GR<-Pop_mean')) {
    met_wide <- all_combi_paths(data_forMA = met_wide, DD = DD, simpleSEM = simpleSEM, Trait = Trait)

    prop_data <- prop_path(data = met_wide, data_MA = data_MA, DD = DD, simpleSEM = simpleSEM)
  }

  trans_allEfS <- met_wide %>%
    tidyr::gather(key, value, -c(Species:ID)) %>%
    tidyr::separate(., key, into = c('Relation', 'Metric'), sep = "/") %>%
    tidyr::spread(., Metric, value)


  ## think if I can add a warning if the specified path is impossible for a given model
  ##  (e.g. Tot_DemRate<-det_Clim for simpleSEM = TRUE)
  subs_merge <- droplevels(data_MA %>%
                             dplyr::distinct(., ID, Country, Continent,
                                             Longitude, Latitude, Taxon,
                                             BirdType, Trait_Categ,
                                             Trait, Demog_rate_Categ,
                                             Demog_rate, Count,
                                             Nyears, WinDur, deltaAIC,
                                             .keep_all = T) %>%
                             subset(.,
                                    select = c(ID, Study_Authors,
                                               Country, Continent,
                                               Longitude, Latitude, Taxon,
                                               BirdType, Trait_Categ, Trait,
                                               Demog_rate_Categ, Demog_rate,
                                               Count, Nyears, WinDur,
                                               deltaAIC, Pvalue, WeathQ,
                                               Ref.day, Ref.month, WindowClose,
                                               LM_std_estimate, LM_std_std.error,
                                               Trait_ageClass, Trend,
                                               GenLength_y_IUCN)))

  tot <- merge(trans_allEfS, subs_merge, by = c('ID'))


  ## subset a specified effect size only
  subs_data <- subset(tot, Relation == Type_EfS)

  # and now get the R2 per relation, to be used on the plot
  trans_allEfS$Response <- unlist(lapply(1:nrow(trans_allEfS), FUN = function(x){
    strsplit(x = trans_allEfS$Relation[x], split = '<')[[1]][1]
  }))

  subs_merge_R2 <- droplevels(data_MA %>%
                                dplyr::distinct(., ID, Country, Continent,
                                                Longitude, Latitude, Taxon,
                                                BirdType, Trait_Categ, Trait,
                                                Demog_rate_Categ, Demog_rate,  Response,
                                                Count, Nyears, WinDur, deltaAIC,
                                                .keep_all = T) %>%
                                subset(.,
                                       select = c(ID, Country, Continent,
                                                  Longitude, Latitude, Taxon,
                                                  BirdType, Trait_Categ, Trait,
                                                  Demog_rate_Categ, Demog_rate,
                                                  Count, Nyears,  Response,
                                                  WinDur, deltaAIC, Pvalue, WeathQ, R.squared,
                                                  LM_std_estimate, LM_std_std.error,
                                                  Trend, Trait_ageClass)))

  tot_R2 <- merge(trans_allEfS, subs_merge_R2, by = c('ID', 'Response'), all.x = TRUE)


  ## subset a specified effect size only
  subs_dataR2 <- subset(tot_R2, Relation == Type_EfS)


  ##  preparing a formula depending on the covariates included
  if(! is.null(Cov_fact)){
    if(! is.null(COV)){
      formul <- paste0('Estimate ~ ', Cov_fact, ' + ', COV, ' - 1')
    } else {
      formul <- paste0('Estimate ~ ', Cov_fact, ' - 1')
    }
  } else {
    if(! is.null(COV)){
      formul <- paste0('Estimate ~ ', COV, ' + 1')
    } else {
      formul <- paste('Estimate ~ 1')
    }
  }

  ## drop the data rows with missing covariate
    if(! is.null(COV)){
    if(length(grep('[+]', COV)) > 0){
      for (j in 1:length(strsplit(COV, '[+]')[[1]])){
        subs_data <- subs_data[! is.na(subs_data[, trimws(strsplit(COV, '[+]')[[1]][j], which = 'both')]), ]
      }
    } else {
      subs_data <- subs_data[! is.na(subs_data[, COV]), ]
      }
    }
    if(! is.null(Cov_fact)){
      subs_data <- subs_data[! is.na(subs_data[, Cov_fact]), ]
    }



  ## fitting the models with ML and REML
  tt.error.ML <- tryCatch(mod_ML <- metafor::rma.mv(stats::as.formula(formul), V = SError^2,
                                                    random = list(~ 1|Species, ~1|ID, ~1|Location),
                                                    data = subs_data,
                                                    method = 'ML', control = list(optimizer = optimize)),
                          error=function(e) e)
  if(is(tt.error.ML,"error")){
    warning(cat('trait category is ', unique(subs_data$Trait_Categ), '\n',
                'demographic rate category is ', unique(subs_data$Demog_rate_Categ), '\n',
                'fitted relation is ', unique(subs_data$Relation), '\n',
                'error when fitting the model with ML \n',
                tt.error.ML[1]$message, '\n'))
  }

  tt.error.REML <- tryCatch(mod_REML <- metafor::rma.mv(stats::as.formula(formul), V = SError^2,
                                                        random = list(~ 1|Species, ~1|ID, ~1|Location),
                                                        data = subs_data,
                                                        method = 'REML',
                                                        control = list(optimizer = optimize)),
                            error=function(e) e)
  if(is(tt.error.REML,"error")){
    warning(cat('trait category is ', unique(subs_data$Trait_Categ), '\n',
                'demographic rate category is ', unique(subs_data$Demog_rate_Categ), '\n',
                'fitted relation is ', unique(subs_data$Relation), '\n',
                'error when fitting the model with REML \n',
                tt.error.REML[1]$message, '\n'))
  }


  ## getting the output depending on whether Cov_fact is included or not (slightly different type of output then
  ## because for covariate we get the estimates per each level but the p value for the whole factor)
  if(! is.null(Cov_fact)){
    if(! is(tt.error.ML,'error') & ! is(tt.error.REML, 'error')){
      out_dat <- data.frame(Levels_Covar = rownames(mod_REML$beta),
                            Estimate = as.numeric(mod_REML$beta),
                            SError = mod_REML$se,
                            EfS_Low = mod_REML$ci.lb,
                            EfS_Upper = mod_REML$ci.ub,
                            AIC_EfS_Covar = rep(AIC(mod_ML), length(as.numeric(mod_REML$beta))),
                            Species.SD = rep(mod_REML$sigma2[grep('Species', mod_REML$s.names)],
                                             length(as.numeric(mod_REML$beta))),
                            ID.SD = rep(mod_REML$sigma2[grep('ID', mod_REML$s.names)],
                                        length(as.numeric(mod_REML$beta))),
                            Location.SD = rep(mod_REML$sigma2[grep('Location', mod_REML$s.names)],
                                              length(as.numeric(mod_REML$beta))))
      if(! is.null(COV)){
        for(i in grep(Cov_fact, rownames(mod_ML$beta))){
          out_dat$z[i] <- stats::anova(mod_ML, btt = grep(Cov_fact, rownames(mod_ML$beta)))$QM
          out_dat$pval_Covar[i] <- stats::anova(mod_ML, btt = grep(Cov_fact, rownames(mod_ML$beta)))$QMp
        }
        if(length(grep('[+]', COV)) > 0){  ## in case COV consists of several elements
          for(j in (length(grep(Cov_fact, rownames(mod_ML$beta))) + 1):
              (length(grep(Cov_fact, rownames(mod_ML$beta))) + length(strsplit(COV, '[+]')[[1]]))){
            out_dat$z[j] <- stats::anova(mod_ML, btt = j)$QM
            out_dat$pval_Covar[j] <- stats::anova(mod_ML, btt = j)$QMp
          }
        } else {
          out_dat$z[nrow(out_dat)] <- stats::anova(mod_ML, btt = nrow(out_dat))$QM
          out_dat$pval_Covar[nrow(out_dat)] <- stats::anova(mod_ML, btt = nrow(out_dat))$QMp
        }
      } else {
        out_dat$z <- rep(stats::anova(mod_ML, btt = grep(Cov_fact, rownames(mod_ML$beta)))$QM,
                         nrow(out_dat))
        out_dat$pval_Covar <- rep(stats::anova(mod_ML, btt = grep(Cov_fact, rownames(mod_ML$beta)))$QMp,
                                  nrow(out_dat))
      }


    }
  } else {
    if(! is(tt.error.ML,'error') & ! is(tt.error.REML, 'error')){
      out_dat <- data.frame(Variable = rownames(mod_REML$beta),
                            Estimate = as.numeric(mod_REML$beta), SError = mod_REML$se,
                            EfS_Low = mod_REML$ci.lb, EfS_Upper = mod_REML$ci.ub,
                            AIC_EfS_across = rep(AIC(mod_ML), length(as.numeric(mod_REML$beta))),
                            Species.SD = mod_REML$sigma2[grep('Species', mod_REML$s.names)],
                            ID.SD = mod_REML$sigma2[grep('ID', mod_REML$s.names)],
                            Location.SD = mod_REML$sigma2[grep('Location', mod_REML$s.names)])

      out_dat$z <- unlist(lapply(1:nrow(out_dat), function(x){
        stats::anova(mod_ML, btt = which(rownames(mod_REML$beta) %in%
                                           rownames(mod_REML$beta)[x]))$QM
      }))
      out_dat$pval_across <- unlist(lapply(1:nrow(out_dat), function(x){
        stats::anova(mod_ML, btt = which(rownames(mod_REML$beta) %in%
                                           rownames(mod_REML$beta)[x]))$QMp
      }))
    }
  }

  ## if both models fitted well (no errors), get the output tibble together
  if(! is(tt.error.ML,'error') & ! is(tt.error.REML, 'error')){
    if(Type_EfS %in% c('Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim',
                       'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim',
                       'Ind_GR<-Pop_mean', 'Tot_GR<-Pop_mean')) {
      out_tib <- tibble::tibble(data = list(out_dat),
                                data_EfS = list(subs_data),
                                data_R2 = list(subs_dataR2),
                                prop_data = list(prop_data))
    } else {
      out_tib <- tibble::tibble(data = list(out_dat),
                                data_EfS = list(subs_data),
                                data_R2 = list(subs_dataR2))
    }
  }

  return(out_tib)
}
