#' Save a table with statistics for correlations among the climate effects on traits
#'
#' \code{tab_CorClimEff} saves the table with the statistics extracted
#' from the model fitted to assess the effect of specific predictors on
#' the correlation between the temperature and precipitation effects on traits.
#'
#' @param data Data frame including for each study the correlation between
#' the effects of precipitaion and temperature on traits, and explanatory variables.
#' The correlation has to be coded as binary (1/0), where 1 = positive, and 0 =
#' negative correlation.
#' @param table_name A string specifying the path and the name of the
#' table to be created.
#' @param explanatory Character indicating the name of the variable whose
#' impact on correlation between the effects of precipitation and temperature on
#' traits is to be tested.
#'
#' @export
#'
#' @return A list containing two objects: (1) a dataframe with the statistics for
#' the effect of a specified predictor on the correlation between the effects of precipitation
#' and temperature on traits; and (2) a model object for the full model. The returned table with
#' statistics is saved (invisibly) as an .xlsx sheet.
#'
#' @examples
#' Coef_ClTr_TP <- readRDS(file = './data/Corr_ClTrEffects.rds')
#' lat_m <- tab_CorClimEff(data = Coef_ClTr_TP,
#'                         table_name = './tables/Latitude_Eff_CorClimEffects',
#'                         explanatory = 'Latitude')
#'
tab_CorClimEff <- function(data, table_name,
                           explanatory){
  formul_full <- paste0('Cor_bin ~ ', explanatory, ' + Trait_Categ +
                        Demog_rate_Categ + Trait_Categ*Demog_rate_Categ + ',
                        explanatory, '*Trait_Categ + ', explanatory,
                        '*Demog_rate_Categ')

  mod_full <- glm(stats::as.formula(formul_full), data = data, family = 'binomial'(link = 'logit'))
  mod_nointDR <- update(mod_full, stats::as.formula(paste(".~.-", explanatory, ':Demog_rate_Categ')))
  mod_nointTr <- update(mod_full, stats::as.formula(paste(".~.-", explanatory, ':Trait_Categ')))
  mod_nointTrDR <- update(mod_full, stats::as.formula(paste(".~.-", 'Trait_Categ:Demog_rate_Categ')))
  mod_noint <- update(mod_nointTrDR, stats::as.formula(paste(".~.-", explanatory, ':Trait_Categ -',
                                                             explanatory, ':Demog_rate_Categ')))
  mod_noDR <- update(mod_noint, stats::as.formula(paste(".~.-", 'Demog_rate_Categ')))
  mod_noTr <- update(mod_noint, stats::as.formula(paste(".~.-", 'Trait_Categ')))
  mod_noExp <- update(mod_noint, stats::as.formula(paste(".~.-", explanatory)))
  mod_incept <- glm(stats::as.formula(paste0('Cor_bin  ~ 1')), data = data, family = 'binomial'(link = 'logit'))
  mod_null <- glm(stats::as.formula(paste0('Cor_bin ~ 0')), data = data, family = 'binomial'(link = 'logit'))


  test_intDR <- anova(mod_nointDR, mod_full, test = 'LRT')
  test_intTr <- anova(mod_nointTr, mod_full, test = 'LRT')
  test_intTrDR <- anova(mod_nointTrDR, mod_full, test = 'LRT')
  test_DR <- anova(mod_noDR, mod_noint, test = 'LRT')
  test_Tr <- anova(mod_noTr, mod_noint, test = 'LRT')
  test_Exp <- anova(mod_noExp, mod_noint, test = 'LRT')
  test_incept <- anova(mod_null, mod_incept, test = 'LRT')

  out <- data.frame(Parameter = c(paste0(explanatory, ':Demog_rate_Categ'),
                                  paste0(explanatory, ':Trait_Categ'),
                                  'Trait_Categ:Demog_rate_Categ', 'Demog_rate_Categ',
                                  'Trait_Categ', paste0(explanatory), 'Intercept'),
                    DF = c(test_intDR$Df[2], test_intTr$Df[2],
                           test_intTrDR$Df[2], test_DR$Df[2], test_Tr$Df[2],
                           test_Exp$Df[2], test_incept$Df[2]),
                    Deviance = c(test_intDR$Deviance[2], test_intTr$Deviance[2],
                              test_intTrDR$Deviance[2], test_DR$Deviance[2],
                              test_Tr$Deviance[2], test_Exp$Deviance[2],
                              test_incept$Deviance[2]),
                    pval = c(test_intDR$`Pr(>Chi)`[2], test_intTr$`Pr(>Chi)`[2],
                             test_intTrDR$`Pr(>Chi)`[2], test_DR$`Pr(>Chi)`[2],
                             test_Tr$`Pr(>Chi)`[2], test_Exp$`Pr(>Chi)`[2],
                             test_incept$`Pr(>Chi)`[2]))

  out$`p value` <- numeric(length = nrow(out))
  for(i in 1:nrow(out)){
    if(! is.na(out$pval[i])){
      if (out$pval[i] < 0.0001){
        out$`p value`[i] <- '<0.0001'
      } else {
        out$`p value`[i] <- format(round(out$pval[i], 4), nsmall = 4, scientific = FALSE)}
    }
  }
  out$pval <- NULL
  out$Deviance <- format(round(out$Deviance, 3), nsmall = 3, scientific = FALSE)
  save_xlsx(table = out,
            table_name = table_name)

  return(list(out, mod_full))
}
