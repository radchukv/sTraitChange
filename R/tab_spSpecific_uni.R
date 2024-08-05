#' Save a table with statistics for species-specific analyses
#'
#' \code{tab_spSpecific} saves the table with the statistics exracted
#' from the model fitted to test the hypothesis about the effect of
#' a specified species characteristic on the indirect path.
#'
#' @param mod_mv Model object of class rma.mv fitted to test the
#' hypothesis about the effect of the specified species characteristic
#' on the indirect path.
#' @param table_name A string specifying the path and the name of the
#' table to be created.
#' @param explanatory Character indicating the name of the variable
#' corresponding to a species characteristic, for which the
#' hypothesis is being tested.
#'
#' @inheritParams fit_mod
#' @export
#'
#' @return A dataframe with the statistics for the test of the effect of
#' a specified characteristic on the indirect path, and saves this dataframe
#' (invisibly) as an .xlsx sheet.
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_overall_sig/PathCoefs_allMods_Temp_AlsoEstimatedRelations.RDS')
#' traits <- read.csv('./data-raw/Species_traits_Subset_11_09.csv')
#' Coef_Aut_IndDR <- droplevels(Coefs_Aut %>%
#' dplyr::filter(., Relation == 'Ind_DemRate<-det_Clim'))
#' Coefs_Aut_sp <- base::merge(Coef_Aut_IndDR, traits_proc, by = 'Species')
#' mod_genLength <- metafor::rma.mv(Estimate ~ GenLength_y_IUCN + Trait_Categ + Demog_rate_Categ +
#'                             GenLength_y_IUCN * Trait_Categ + GenLength_y_IUCN * Demog_rate_Categ +
#'                             Pvalue,
#'                             V = SError^2, random = list(~ 1|Species, ~1|ID, ~1|Location),
#'                             data = Coefs_Aut_sp, method = 'ML')
#' tab_spSpecific(mod_mv = mod_genLength, table_name = './tables/GenLength_Temp', explanatory = 'GenLength_y_IUCN')
#'
tab_spSpecific_uni <- function(mod_mv, table_name,
                           explanatory, interact_fac){
  stats <- stats::coef(summary(mod_mv))
  stats$Parameter <- rownames(stats)
  stats <- stats %>%
    dplyr::select(., -c(ci.lb, ci.ub))
  colnames(stats) <- c('Estimate', 'SE', 'Chi2', 'pval', 'Parameter')

  ## estimate the effect for an interaction
  Inter <- stats::anova(mod_mv, btt = grep(interact_fac, stats$Parameter))

  stats$DF <- rep(1, nrow(stats))
  ## replace the statistics with the  performed omnibus tests
  stats <- replace_stats(data = stats, variable = interact_fac, stats_out = Inter)

  ## 4. for explanatory if it is not continuous - actually the stats on Chi2 has to be updated anyways
  for(i in 1:length(explanatory)){
    Explan <- stats::anova(mod_mv, btt = which(stats$Parameter %in%
                                              stats$Parameter[grepl(explanatory[i], stats$Parameter) &
                                                                !grepl(':', stats$Parameter)]))
    stats <- replace_stats(data = stats, variable = explanatory[i], stats_out = Explan)
  }


  stats <- stats[c('Parameter', 'Estimate', 'SE', 'Chi2', 'pval', 'DF')]
  stats$`p value` <- numeric(length = nrow(stats))
  for(i in 1:nrow(stats)){
    if(! is.na(stats$pval[i])){
      if (stats$pval[i] < 0.0001){
        stats$`p value`[i] <- '<0.0001'
      } else {
        stats$`p value`[i] <- format(round(stats$pval[i], 4), nsmall = 4, scientific = FALSE)}
    }
  }
  stats$pval <- NULL
  stats$Estimate <- format(round(stats$Estimate, 3), nsmall = 3, scientific = FALSE)
  stats$SE <- format(round(stats$SE, 3), nsmall = 3, scientific = FALSE)
  stats$Chi2 <- format(round(stats$Chi2, 3), nsmall = 3, scientific = FALSE)
  stats$`p value`[stats$`p value` == 0] <- NA
  save_xlsx(table = stats,
            table_name = table_name)

  return(stats)
}
