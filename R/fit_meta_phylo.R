#' Fit a meta-analytical model to assess the global effects across the studies
#'
#' \code{fit_meta_phylo} fits mixed-effects meta-analytical model to extract a global
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
#' @param A A variance-covariance matrix based on phylogeny.
#' @param des.matrix A character specifying what type of design matrix to use if the factor
#' is included as one of the explanatory variables. By default treatment contrasts is used
#'  (des.matrix = "treatm.contrasts") but the user can also request identity matrix by
#'  specifying "identity".
#'
#' @inheritParams fit_mod
#' @importFrom magrittr "%>%"
#'
#' @export
#'
#' @return Returns a tibble that includes five data frames:
#' 1. 'data' - the data frame with the estimate
#' of the global effect size, its standard error, its significance, the AIC of the fitted
#' mixed-effects model and the variances explained by random effects;
#' 2. 'data_Efs' - the data frame with the effect sizes and their standard errors for each study;
#' 3. 'heter_mod' - the data frame with estimates of heterogeneity
#' typically extracted for meta-analysis;
#' 4. ML_mod - results of the meta-analysis correcting for phylogeny (if that
#' improves mode lit to the data) and fitted with ML;
#' 5. REML_mod - results of the meta-analysis correcting for phylogeny (if that
#' improves mode lit to the data) and fitted with REML.
#' The AIC and the p values are obtained from the models fitted with ML, whereas
#' the parameter estimates from the models fitted using REML.
#'
#' @examples
#' # prepare dataset, select only studies with phenological traits
#' Coefs_phenClim <- subset(dataPaths, Relation == 'Trait_mean<-det_Clim' &
#' Trait_Categ == 'Phenological')
#' Coefs_phenClim <- Coefs_phenClim %>%
#'                   dplyr::mutate(Species = dplyr::case_when(
#'                          Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
#'                          Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
#'                          Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
#'                          Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
#'                          Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
#'                          Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
#'                          TRUE ~ Species))
#'
#' Coefs_phenClim$Species <- unlist(lapply(1:nrow(Coefs_phenClim), FUN = function(x){
#'   binary <- strsplit(as.character(Coefs_phenClim$Species[x]), " ")
#'   Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")}))
#' Coefs_phenClim$Sp_phylo <- Coefs_phenClim$Species
#'
#' test_noCovar <- fit_meta_phylo(data_MA = Coefs_phenClim,
#'                               Type_EfS = 'Trait_mean<-det_Clim',
#'                               Cov_fact = NULL, COV = NULL,
#'                               DD = 'n_effectGR',
#'                               simpleSEM = TRUE,
#'                               A = phyloMat)
#' # meta-analysis with a quantitative covariate
#' test_PValCovar <- fit_meta_phylo(data_MA = Coefs_phenClim,
#'                               Type_EfS = 'Trait_mean<-det_Clim',
#'                               Cov_fact = NULL, COV = 'Pvalue',
#'                               DD = 'n_effectGR',
#'                               simpleSEM = TRUE,
#'                               A = phyloMat)
#'
#' test_WeathQ <- fit_meta_phylo(data_MA = Coefs_phenClim,
#'                                    Type_EfS = 'Trait_mean<-det_Clim',
#'                                    Cov_fact = 'WeathQ',
#'                                    COV = NULL,
#'                                    DD = 'n_effectGR',
#'                                    simpleSEM = TRUE,
#'                                    A = phyloMat)
#' test_WeathQ
fit_meta_phylo <- function(data_MA, Type_EfS = 'Trait_mean<-det_Clim',
                     Cov_fact = NULL, COV = NULL, optimize = 'uobyqa',
                     DD = 'n_effectDGR', simpleSEM = FALSE,
                     Trait = FALSE, A = Mat_phylo, des.matrix = "treatm.contrasts"){
  ## calculating indirect effects, total effects and their SEs
  forTrans <- subset(data_MA, select = c(Estimate,  Std.Error, Relation, Species, Sp_phylo, Location, ID))
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
  }

  trans_allEfS <- met_wide %>%
    tidyr::gather(key, value, -c(Species:ID)) %>%
    tidyr::separate(., key, into = c('Relation', 'Metric'), sep = "/") %>%
    tidyr::spread(., Metric, value)


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
                                               Trait_ageClass,
                                               GenLength_y_IUCN)))

  tot <- merge(trans_allEfS, subs_merge, by = c('ID'))


  ## subset a specified effect size only
  subs_data <- subset(tot, Relation == Type_EfS)

  ##  preparing a formula depending on the covariates included
  if(des.matrix == 'treatm.contrasts') {
  if(! is.null(Cov_fact)){
    if(! is.null(COV)){
      formul <- paste0('Estimate ~ ', Cov_fact, ' + ', COV)
    } else {
      formul <- paste0('Estimate ~ ', Cov_fact)
    }
  } else {
    if(! is.null(COV)){
      formul <- paste0('Estimate ~ ', COV, ' + 1')
    } else {
      formul <- paste('Estimate ~ 1')
    }
  }
  } else if (des.matrix == 'identity') {
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

  ## fitting the models without phylogeny with ML and REML
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


  ## fitting the models with phylogeny (and non-phylogeny related effect of species)
  ## with ML and REML
  tt.error.phylo.ML <- tryCatch(mod_phylo_ML <- metafor::rma.mv(stats::as.formula(formul), V = SError^2,
                                                                random = list(~ 1|Species, ~1|ID, ~1|Location,
                                                                              ~1|Sp_phylo),
                                                                R = list(Sp_phylo = A),
                                                                data = subs_data,
                                                                method = 'ML', control = list(optimizer = optimize)),
                                error=function(e) e)
  if(is(tt.error.phylo.ML,"error")){
    warning(cat('trait category is ', unique(subs_data$Trait_Categ), '\n',
                'demographic rate category is ', unique(subs_data$Demog_rate_Categ), '\n',
                'fitted relation is ', unique(subs_data$Relation), '\n',
                'error when fitting the model with ML and phylogeny \n',
                tt.error.phylo.ML[1]$message, '\n'))
  }

  tt.error.phylo.REML <- tryCatch(mod_phylo_REML <- metafor::rma.mv(stats::as.formula(formul), V = SError^2,
                                                                    random = list(~ 1|Species, ~1|ID, ~1|Location,
                                                                                  ~1|Sp_phylo),
                                                                    R = list(Sp_phylo = A),
                                                                    data = subs_data,
                                                                    method = 'REML',
                                                                    control = list(optimizer = optimize)),
                                  error=function(e) e)
  if(is(tt.error.phylo.REML,"error")){
    warning(cat('trait category is ', unique(subs_data$Trait_Categ), '\n',
                'demographic rate category is ', unique(subs_data$Demog_rate_Categ), '\n',
                'fitted relation is ', unique(subs_data$Relation), '\n',
                'error when fitting the model with REML and phylogeny \n',
                tt.error.phylo.REML[1]$message, '\n'))
  }

  if(AIC(mod_REML) > AIC(mod_phylo_REML)){
    sel_mod_REML <- mod_phylo_REML
    sel_mod_ML <- mod_phylo_ML
    tt.error.sel.ML <- tt.error.phylo.ML
    tt.error.sel.REML <- tt.error.phylo.REML
  } else {
    sel_mod_REML <- mod_REML
    sel_mod_ML <- mod_ML
    tt.error.sel.ML <- tt.error.ML
    tt.error.sel.REML <- tt.error.REML
  }
  ## getting the output depending on whether Cov_fact is included or not (slightly different type of output then
  ## because for covariate we get the estimates per each level but the p value for the whole factor)
  if(! is.null(Cov_fact)){
    if(! is(tt.error.sel.ML,'error') & ! is(tt.error.sel.REML, 'error')){
      out_dat <- data.frame(Levels_Covar = rownames(sel_mod_REML$beta),
                            Estimate = as.numeric(sel_mod_REML$beta),
                            SError = sel_mod_REML$se,
                            EfS_Low = sel_mod_REML$ci.lb,
                            EfS_Upper = sel_mod_REML$ci.ub,
                            AIC_EfS_Covar = rep(AIC(sel_mod_ML), length(as.numeric(sel_mod_REML$beta))),
                            AIC_phylo = rep(AIC(mod_phylo_REML), length(as.numeric(sel_mod_REML$beta))),
                            AIC_nophyl = rep(AIC(mod_REML), length(as.numeric(sel_mod_REML$beta))),
                            Species.SD = rep(sel_mod_REML$sigma2[grep('Species', sel_mod_REML$s.names)],
                                             length(as.numeric(sel_mod_REML$beta))),
                            ID.SD = rep(sel_mod_REML$sigma2[grep('ID', sel_mod_REML$s.names)],
                                        length(as.numeric(sel_mod_REML$beta))),
                            Location.SD = rep(sel_mod_REML$sigma2[grep('Location', sel_mod_REML$s.names)],
                                              length(as.numeric(sel_mod_REML$beta))),
                            Phylo.SD = rep(mod_phylo_REML$sigma2[grep('Sp_phylo', mod_phylo_REML$s.names)],
                                             length(as.numeric(sel_mod_REML$beta))))

      if(des.matrix == "treatm.contrasts") {
        out_dat$Chi2[1] <- stats::anova(sel_mod_ML, btt = 1)$QM
        out_dat$pval_Covar[1] <- stats::anova(sel_mod_ML, btt = 1)$QMp
        for(i in grep(Cov_fact, rownames(sel_mod_ML$beta))){
          out_dat$Chi2[i] <- stats::anova(sel_mod_ML, btt = c(1, grep(Cov_fact, rownames(sel_mod_ML$beta))))$QM
          out_dat$pval_Covar[i] <- stats::anova(sel_mod_ML, btt = c(1, grep(Cov_fact, rownames(sel_mod_ML$beta))))$QMp
        }
      } else if(des.matrix == "identity") {
        for(i in grep(Cov_fact, rownames(sel_mod_ML$beta))){
          out_dat$Chi2[i] <- stats::anova(sel_mod_ML, btt = grep(Cov_fact, rownames(sel_mod_ML$beta)))$QM
          out_dat$pval_Covar[i] <- stats::anova(sel_mod_ML, btt = grep(Cov_fact, rownames(sel_mod_ML$beta)))$QMp
        }
      }

      if(! is.null(COV)){
        if(length(grep('[+]', COV)) > 0){  ## in case COV consists of several elements
          for(j in (length(grep(Cov_fact, rownames(sel_mod_ML$beta))) + 1):
              (length(grep(Cov_fact, rownames(sel_mod_ML$beta))) + length(strsplit(COV, '[+]')[[1]]))){
            out_dat$Chi2[j] <- stats::anova(sel_mod_ML, btt = j)$QM
            out_dat$pval_Covar[j] <- stats::anova(sel_mod_ML, btt = j)$QMp
          }
        } else {
          out_dat$Chi2[nrow(out_dat)] <- stats::anova(sel_mod_ML, btt = nrow(out_dat))$QM
          out_dat$pval_Covar[nrow(out_dat)] <- stats::anova(sel_mod_ML, btt = nrow(out_dat))$QMp
        }
      }
    }
  } else {  # there is no factor as an explanatory
    if(! is(tt.error.sel.ML,'error') & ! is(tt.error.sel.REML, 'error')){
      out_dat <- data.frame(Variable = rownames(sel_mod_REML$beta),
                            Estimate = as.numeric(sel_mod_REML$beta), SError = sel_mod_REML$se,
                            EfS_Low = sel_mod_REML$ci.lb, EfS_Upper = sel_mod_REML$ci.ub,
                            AIC_EfS_Covar = rep(AIC(sel_mod_ML), length(as.numeric(sel_mod_REML$beta))),
                            AIC_phylo = rep(AIC(mod_phylo_REML), length(as.numeric(sel_mod_REML$beta))),
                            AIC_nophyl = rep(AIC(mod_REML), length(as.numeric(sel_mod_REML$beta))),
                            Species.SD = sel_mod_REML$sigma2[grep('Species', sel_mod_REML$s.names)],
                            ID.SD = sel_mod_REML$sigma2[grep('ID', sel_mod_REML$s.names)],
                            Location.SD = sel_mod_REML$sigma2[grep('Location', sel_mod_REML$s.names)],
                            Phylo.SD = rep(mod_phylo_REML$sigma2[grep('Sp_phylo', mod_phylo_REML$s.names)],
                                           length(as.numeric(sel_mod_REML$beta))))

      out_dat$Chi2 <- unlist(lapply(1:nrow(out_dat), function(x){
        stats::anova(sel_mod_ML, btt = which(rownames(sel_mod_REML$beta) %in%
                                               rownames(sel_mod_REML$beta)[x]))$QM
      }))
      out_dat$pval_Covar <- unlist(lapply(1:nrow(out_dat), function(x){
        stats::anova(sel_mod_ML, btt = which(rownames(sel_mod_REML$beta) %in%
                                               rownames(sel_mod_REML$beta)[x]))$QMp
      }))
    }
  }

  ## if both models fitted well (no errors), get the output tibble together
  if(! is(tt.error.sel.ML,'error') & ! is(tt.error.sel.REML, 'error')){
    het_mod <- get_heterog(mod= tt.error.sel.REML, data = subs_data)
    if(Type_EfS %in% c('Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim',
                       'Tot_DemRate<-det_Clim', 'Tot_GR<-det_Clim',
                       'Ind_GR<-Pop_mean', 'Tot_GR<-Pop_mean')) {
      out_tib <- tibble::tibble(data = list(out_dat),
                                data_EfS = list(subs_data),
                                heter_mod = list(het_mod),
                                ML_mod = list(tt.error.sel.ML),
                                REML_mod = list(tt.error.sel.REML))
    } else {
      out_tib <- tibble::tibble(data = list(out_dat),
                                data_EfS = list(subs_data),
                                heter_mod = list(het_mod),
                                ML_mod = list(tt.error.sel.ML),
                                REML_mod = list(tt.error.sel.REML))
    }
  }

  return(out_tib)
}
