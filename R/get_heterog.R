#' Assess model heterogeneity
#'
#' \code{get_heterog} computes the heterogeneity metrics for meta-analysis.
#'
#' This function computes the heterogeneity metrics for the meta-analysis,
#' including I2, Q and p value associated with !
#'
#' @param mod A fitted mixed-effects meta-analytical model.
#' @param data Data frame used to fit the meta-analytical model,
#' which contains effect sizes and the respective sampling variances.
#'
#' @export
#' @return A tibble of the calculated heterogeneity metrics: Q reflects
#' the total amount of heterogeneity, Qp is the p value showing whether the observed
#' heterogeneity is larger than would be expected based on sampling variability,
#' I2 reflects the proportion of the total heterogeneity due to between-study
#' variance, ranges from 0 to 1, and I2_perRand is the proportion of the total
#' heterogeneity per random factor included in the model.
#' @examples
#' # prepare the data to fit the model
#' Coefs_phenClim <- subset(dataPaths, Relation == 'Trait_mean<-det_Clim' &
#' Trait_Categ == 'Phenological')
#'  forTrans <- subset(Coefs_phenClim, select = c(Estimate,  Std.Error, Relation,
#'  Species, Location, ID))
#'  forTrans <- forTrans %>%
#'              dplyr::rename(SError = Std.Error)
#'  met_wide <- forTrans %>%
#'      tidyr::gather(variable, value, -(Relation:ID)) %>%
#'      tidyr::unite(temp, Relation, variable, sep = '/') %>%
#'      tidyr::spread(temp, value)
#' trans_allEfS <- met_wide %>%
#'     tidyr::gather(key, value, -c(Species:ID)) %>%
#'     tidyr::separate(., key, into = c('Relation', 'Metric'), sep = "/") %>%
#'     tidyr::spread(., Metric, value)
#' subs_merge <- droplevels(Coefs_phenClim %>%
#'         dplyr::distinct(., ID, Country, Continent,
#'         Longitude, Latitude, Taxon,
#'         BirdType, Trait_Categ,
#'         Trait, Demog_rate_Categ,
#'         Demog_rate, Count,
#'         Nyears, WinDur, deltaAIC,
#'         .keep_all = TRUE) %>%
#'         subset(., select = c(ID, Study_Authors,
#'         Country, Continent,
#'         Longitude, Latitude, Taxon,
#'         BirdType, Trait_Categ, Trait,
#'         Demog_rate_Categ, Demog_rate,
#'         Count, Nyears, WinDur,
#'         deltaAIC, Pvalue, WeathQ,
#'         Ref.day, Ref.month, WindowClose,
#'         Trait_ageClass, GenLength_y_IUCN)))
#' tot <- merge(trans_allEfS, subs_merge, by = c('ID'))
#' # subset a specified effect size only
#' subs_data <- subset(tot, Relation == 'Trait_mean<-det_Clim')
#'
#'# fit model
#' mod_REML <- metafor::rma.mv(Estimate ~ 1, V = SError^2,
#'                             random = list(~ 1|Species, ~1|ID,
#'                             ~1|Location),
#'                             data = subs_data,
#'                             method = 'REML')
#' het_mod <- get_heterog(mod= mod_REML, data = subs_data)
#'
get_heterog <- function(mod, data){
  W <- diag(1/data$SError^2)
  X <- stats::model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 <- sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))
  I2_perRand <- mod$sigma2 / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))
  return(tibble::tibble(Q = mod$QE, Qp = mod$QEp,
                        I2 = I2, I2perRand = I2_perRand))
}
