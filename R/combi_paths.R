#' Provide estimates of median and SE of an indirect (focal) pathway by bootstrap
#'
#' \code{ind_path} estimates median and SE of an indirect pathway
#' by bootstrap
#'
#' @param x Numeric reflecting the first path coefficient estimate
#' in the indirect path
#' @param y Numeric reflecting the second path coefficient estimate
#' in the indirect path
#' @param z Numeric reflecting the third path coefficient estimate
#'  in the indirect path, defaults to NULL, in which case
#'  the indirect path consists of two pathways.
#' @param x.se Numeric reflecting the standard error for
#'  the first path coefficient estimate.
#' @param y.se Numeric reflecting the standard error for
#' the second path coefficient estimate.
#' @param z.se Numeric reflecting the standard error for
#' the third path coefficient estimate, defaults to NULL for
#' the case of the indirect path consisting of two pathways.
#' @param numrep Numeric for the number of replicates to use in a
#' bootstrap.
#' @param omega Numeric reflecting the fourth path coefficient estimate
#' (applicable only if Trait = TRUE in \code{fit_SEM}). Defaults to NULL.
#' @param omega.se Numeric reflecting the standard errir for the
#' fourth path coefficient estimate (applicable only if T
#' rait = TRUE in \code{fit_SEM}). Defaults to NULL.
#'
#' @return A dataframe that contains four columns: a median of the
#' bootstrap values as an estimate of an indirect path coefficient,
#' lower and upper confidence intervals, and SE of the estimate.
#' @export
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_forSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#'   forTrans <- subset(Coefs_Aut,
#'               select = c(Estimate,  Std.Error, Relation,
#'               Species, Location, ID)) %>%
#'               dplyr::rename(SError = Std.Error)
#' met_wide <- forTrans %>%
#'   tidyr::gather(variable, value, -(Relation:ID)) %>%
#'   tidyr::unite(temp, Relation, variable, sep = '/') %>%
#'   tidyr::spread(temp, value)
#'
#' ind_p_Dem <- ind_path(x = met_wide$`Trait_mean<-det_Clim/Estimate`[1],
#'                   y = met_wide$`Demog_rate_mean<-Trait_mean/Estimate`[1],
#'                   z = NULL,
#'                   x.se = met_wide$`Trait_mean<-det_Clim/SError`[1],
#'                   y.se = met_wide$`Demog_rate_mean<-Trait_mean/SError`[1],
#'                   z.se = NULL, numrep = 10000)
#' ind_p_GR <- ind_path(x = met_wide$`Trait_mean<-det_Clim/Estimate`[2],
#'                      y = met_wide$`Demog_rate_mean<-Trait_mean/Estimate`[2],
#'                      z = met_wide$`GR<-Demog_rate_mean/Estimate`[2],
#'                      x.se = met_wide$`Trait_mean<-det_Clim/SError`[2],
#'                      y.se = met_wide$`Demog_rate_mean<-Trait_mean/SError`[2],
#'                      z.se = met_wide$`GR<-Demog_rate_mean/SError`[2],
#'                      numrep = 10000)
#'
ind_path <- function(x, y, z = NULL,
                     x.se, y.se,
                     z.se = NULL,
                     numrep = 10000,
                     omega = NULL,
                     omega.se = NULL){
  if(! is.null(z)){
    path = replicate(numrep, stats::rnorm(1, x, x.se) * stats::rnorm(1, y, y.se) *
                       stats::rnorm(1, z, z.se) +
                       stats::rnorm(1, x, x.se) * stats::rnorm(1, omega, omega.se))

  } else {
    path = replicate(numrep, stats::rnorm(1, x, x.se) * stats::rnorm(1, y, y.se))
  }
  med <- stats::median(path)
  lCI <- stats::quantile(path, 0.025)
  uCI <- stats::quantile(path, 0.975)
  if ((lCI > 0 & uCI > 0) | (lCI < 0 & uCI < 0)){
    SE = (uCI - lCI) / 4
    } else {
      SE = (abs(uCI) + abs(lCI)) / 4
    }
  return(data.frame(Median = med, lCI = lCI,
                    uCI = uCI, SE = SE))
}


#' Provide estimates of median and SE of the total pathway by bootstrap
#'
#' \code{tot_path} estimates median and SE of the total pathway
#' by bootstrap
#'
#' @param direct Numeric reflecting the estimate of the direct
#' path coefficient.
#' @param indir Numeric reflecting the estimate of the indirect
#' path coefficient.
#' @param ClDem Numeric reflecting the estimate of the path from
#' climate to demographic rate, in case the total path from
#' climate to growth rate is estimated. Defaults to NULL.
#' @param DemGR Numeric reflecting the estimate of the path from
#' demographic rate to growth rate, in case the total path from
#' climate to growth rate is estimated. Defaults to NULL.
#' @param ClDem.se Numeric reflecting the standard error for the path from
#' climate to demographic rate.
#' @param DemGR.se Numeric reflecting the standard error for the path from
#' demographic rate to growth rate.
#' @param direct.se Numeric reflecting the standard error for the direct
#' path coefficient.
#' @param indir.se Numeric reflecting the standard error for the indirect
#' path coefficient.
#' @param numrep Numeric for the number of replicates to use in a
#' bootstrap.
#'
#' @return A dataframe that contains four columns: a median of the
#' bootstrap values as an estimate of an indirect path coefficient,
#' lower and upper confidence intervals, and SE of the estimate.
#' @export
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_forSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#'   forTrans <- subset(Coefs_Aut,
#'               select = c(Estimate,  Std.Error, Relation,
#'               Species, Location, ID)) %>%
#'               dplyr::rename(SError = Std.Error)
#' met_wide <- forTrans %>%
#'   tidyr::gather(variable, value, -(Relation:ID)) %>%
#'   tidyr::unite(temp, Relation, variable, sep = '/') %>%
#'   tidyr::spread(temp, value)
#'
#' ind_p_Dem <- ind_path(x = met_wide$`Trait_mean<-det_Clim/Estimate`[1],
#'                   y = met_wide$`Demog_rate_mean<-Trait_mean/Estimate`[1],
#'                   z = NULL,
#'                   x.se = met_wide$`Trait_mean<-det_Clim/SError`[1],
#'                   y.se = met_wide$`Demog_rate_mean<-Trait_mean/SError`[1],
#'                   z.se = NULL, numrep = 10000)
#' tot_Dem <- tot_path(direct = met_wide$`Demog_rate_mean<-det_Clim/Estimate`[1],
#'                     indir = ind_p_Dem$Median,
#'                     ClDem = NULL,  DemGR= NULL,
#'                     direct.se = met_wide$`Demog_rate_mean<-det_Clim/SError`[1],
#'                     indir.se = ind_p_Dem$SE,
#'                     ClDem.se = NULL, DemGR.se = NULL,
#'                     numrep = 10000)
#'
tot_path <- function(direct, indir,
                     ClDem = NULL,  DemGR= NULL,
                     direct.se, indir.se,
                     ClDem.se = NULL, DemGR.se = NULL,
                     numrep = 10000){
  if(! is.null(ClDem)){
    tot = replicate(numrep, stats::rnorm(1, indir, indir.se) +
                      stats::rnorm(1, direct, direct.se) +
                      stats::rnorm(1, ClDem, ClDem.se) * stats::rnorm(1, DemGR, DemGR.se))
  } else {
    tot = replicate(numrep, stats::rnorm(1, indir, indir.se) +
                      stats::rnorm(1, direct, direct.se))
  }
  med <- stats::median(tot)
  lCI <- stats::quantile(tot, 0.025)
  uCI <- stats::quantile(tot, 0.975)
  if ((lCI > 0 & uCI > 0) | (lCI < 0 & uCI < 0)){
    SE = (uCI - lCI) / 4
  } else {
    SE = (abs(uCI) + abs(lCI)) / 4
  }
  return(data.frame(Median = med, lCI = lCI,
                    uCI = uCI, SE = SE))
}
