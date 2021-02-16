#' Fit beta regression model to the proportional contribution data
#'
#' \code{fit_mod_beta} fits beta regression model either with or
#' without an additional explanatory (species trait), as requested
#'
#' @param prop_data Data frame containing for each study ID proportional
#' contribution of the trait-mediated effect either to the demographic rate
#' or to the growth rate (previously prepared by the user)
#' @param explanatory Character specifying the name of the additional explanatory
#' (i.e. species-specific trait) to be included in the model.
#' @inheritParams save_xlsx
#'
#' @export
#'
#' @return A list of length 2, with the first object containing
#' the fit statistics for each explanatory in the model (performed
#' with Wald test), and the secon dobject containing the fitted full model.
#'
#' @examples
#' prop_T_simple <- readRDS(file = './output_all_simpleSEM/all_prop_EffectTempOnTrait_GroupsTraitDem.RDS')
#' prop_T_simple <- prop_T_simple %>%
#'       dplyr::mutate(., Climatic_var = 'Temperature')
#' prop_P_simple <- readRDS(file = './output_all_simpleSEM/all_prop_EffectPrecipOnTrait_GroupsTraitDem.RDS')
#' prop_P_simple <- prop_P_simple %>%
#'      dplyr::mutate(., Climatic_var = 'Precipitation')
#' prop_sim <- rbind(prop_T_simple, prop_P_simple)
#' prop_long <- prop_sim %>%
#'      tidyr::gather(key = 'Path', value = 'Proportion',
#'      dplyr::contains('prop_'))
#' prop_IndGR <- prop_long %>%
#'      dplyr::filter(., Path == 'prop_ind_GR')
#' propGR_bet <- fit_mod_beta(prop_data = prop_IndGR,
#'                            table_name = './tables/PropGR_general',
#'                            explanatory = NULL)
#' propGR_bet[[1]]
#' mod_GR <- propGR_bet[[2]]
#' mod_GR
#'
fit_mod_beta <- function(prop_data, table_name,
                          explanatory= NULL){
  if(! is.null(explanatory)){
    ## to start with including all the variables (also clim var, which was not signif. in the model without the explanatory)
    formul_full <- paste0('Proportion ~ ', explanatory, ' + Climatic_var + Trait_Categ +
                          Climatic_var:Trait_Categ +', explanatory, ':Climatic_var + ',
                          explanatory, ':Trait_Categ')
  }else {
    formul_full <- paste0('Proportion ~ Climatic_var + Trait_Categ + Climatic_var:Trait_Categ')
  }

  mod_full <- betareg::betareg(stats::as.formula(formul_full), data = prop_data)

  ## models fitted no matter whether explnatory is included or not
  mod_nointClTr <- update(mod_full, stats::as.formula(paste(".~.-", 'Climatic_var:Trait_Categ')))
  test_intClTr <- lmtest::waldtest(mod_nointClTr, mod_full)

  if(! is.null(explanatory)){
    ## models unique in case explanatory is included
    mod_nointExpCl <- update(mod_full, stats::as.formula(paste0(".~.- ", explanatory, ':Climatic_var')))
    mod_nointExpTr <- update(mod_full, stats::as.formula(paste0(".~.- ", explanatory, ':Trait_Categ')))

    mod_noint <- update(mod_nointClTr,
                        stats::as.formula(paste0(".~.- ", explanatory, ':Climatic_var - ',
                                                 explanatory, ':Trait_Categ')))
    mod_noExp <- update(mod_noint, stats::as.formula(paste0(".~.- ", explanatory)))

    test_intExpCl <- lmtest::waldtest(mod_nointExpCl, mod_full)
    test_intExpTr <- lmtest::waldtest(mod_nointExpTr, mod_full)
    test_Exp <- lmtest::waldtest(mod_noExp, mod_noint)

  }
  else {
    mod_noint <- mod_nointClTr
  }

  mod_noTr <- update(mod_noint, stats::as.formula(paste(".~.-", 'Trait_Categ')))
  mod_noCl <- update(mod_noint, stats::as.formula(paste(".~.-", 'Climatic_var')))

  mod_incept <- betareg::betareg(stats::as.formula(paste0('Proportion  ~ 1')), data = prop_data)
  mod_null <- betareg::betareg(stats::as.formula(paste0('Proportion ~ 0')), data = prop_data)


  test_Tr <- lmtest::waldtest(mod_noTr, mod_noint)
  test_Cl <- lmtest::waldtest(mod_noCl, mod_noint)
  test_incept <- lmtest::waldtest(mod_null, mod_incept)

  if(! is.null(explanatory)){
    out = data.frame(Parameter = c('Climatic_var:Trait_Categ',
                                   paste0(explanatory, ':Climatic_var'),
                                   paste0(explanatory, ':Trait_Categ'),
                                   'Trait_Categ', 'Climatic_var',
                                   paste0(explanatory), 'Intercept'),
                     DF = c(test_intClTr$Df[2], test_intExpCl$Df[2],
                            test_intExpTr$Df[2], test_Tr$Df[2],
                            test_Cl$Df[2], test_Exp$Df[2], test_incept$Df[2]),
                     Chi2 = c(test_intClTr$Chisq[2], test_intExpCl$Chisq[2],
                              test_intExpTr$Chisq[2], test_Tr$Chisq[2],
                              test_Cl$Chisq[2], test_Exp$Chisq[2],
                              test_incept$Chisq[2]),
                     pval = c(test_intClTr$`Pr(>Chisq)`[2], test_intExpCl$`Pr(>Chisq)`[2],
                              test_intExpTr$`Pr(>Chisq)`[2], test_Tr$`Pr(>Chisq)`[2],
                              test_Cl$`Pr(>Chisq)`[2], test_Exp$`Pr(>Chisq)`[2],
                              test_incept$`Pr(>Chisq)`[2]))
  } else {
    out <- data.frame(Parameter = c('Climatic_var:Trait_Categ', 'Trait_Categ',
                                    'Climatic_var','Intercept'),
                      DF = c(test_intClTr$Df[2], test_Tr$Df[2],
                             test_Cl$Df[2], test_incept$Df[2]),
                      Chi2 = c(test_intClTr$Chisq[2], test_Tr$Chisq[2],
                               test_Cl$Chisq[2], test_incept$Chisq[2]),
                      pval = c(test_intClTr$`Pr(>Chisq)`[2], test_Tr$`Pr(>Chisq)`[2],
                               test_Cl$`Pr(>Chisq)`[2], test_incept$`Pr(>Chisq)`[2]))
  }
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
  out$Chi2 <- format(round(out$Chi2, 3), nsmall = 3, scientific = FALSE)
  save_xlsx(table = out,
            table_name = table_name)

  return(list(out, mod_full))

}
