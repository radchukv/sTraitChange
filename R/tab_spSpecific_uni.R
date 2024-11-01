#' Save a table with statistics for species-specific analyses
#'
#' \code{tab_spSpecific_uni} saves the table with the statistics extracted
#' from the model fitted to test hypothesis about the effects of
#' specified species characteristic(s) on the path from SEMs.
#'
#' @param mod_mv Model object of class rma.mv fitted to test the
#' hypothesis about the effects of the specified species characteristic(s)
#' on a certain path from SEMs.
#' @param table_name A string specifying the path and the name of the
#' table to be created.
#' @param explanatory A character string indicating the name(s) of
#' the variable(s) corresponding to (a) species characteristic(s),
#' for which the hypothesis is being tested.
#' @param interact_fac A character specifying the name of the factor
#' variable with which interactions of other variable(s) are tested,
#' defaults to NULL.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @return A dataframe with the statistics for the test of the effect(s) of
#' (a) specified characteristic(s) on the specific path from SEMs,
#' and saves this dataframe (invisibly) as an .xlsx sheet.
#'
#' @examples
#' dataPaths_phen <- dataPaths %>%
#'     dplyr::filter(Trait_Categ == 'Phenological') %>%
#'     dplyr::mutate(absLat = abs(Latitude))
#'
#' mod_CZ_PhenT_AbsLat <- metafor::rma.mv(Estimate ~ absLat + Pvalue,
#'                               V = Std.Error^2, random = list(~ 1|Species, ~1|ID, ~1|Location),
#'                               data = dataPaths_phen, method = 'ML')
#' summary(mod_CZ_PhenT_AbsLat)
#' tab_spSpecific_uni(mod_mv = mod_CZ_PhenT_AbsLat,
#'                    table_name = './tables/CZ_PhenT_AbsLat',
#'                    explanatory = c('absLat', 'Pvalue'),
#'                    interact_fac = NULL)
#'
tab_spSpecific_uni <- function(mod_mv, table_name,
                           explanatory, interact_fac = NULL){
  stats <- stats::coef(summary(mod_mv))
  stats$Parameter <- rownames(stats)
  stats$DF <- numeric(length = nrow(stats))
  stats <- stats %>%
    dplyr::select(-c(.data$ci.lb, .data$ci.ub))
  colnames(stats) <- c('Estimate', 'SE', 'Chi2', 'pval', 'Parameter', 'DF')

  if(! is.null(interact_fac)){

  ## estimate the effect for an interaction
  Inter <- stats::anova(mod_mv, btt = grep(interact_fac, stats$Parameter))

  stats$DF <- rep(1, nrow(stats))
  ## replace the statistics with the  performed omnibus tests
  stats <- replace_stats(data = stats, variable = interact_fac, stats_out = Inter)
  }

  ## 4. for explanatory if it is not continuous or a factor with 2 levels
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
