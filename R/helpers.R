#' Impute average of the values before and after the missing one
#'
#' \code{imput_ma} imputes the mean of the values before and after
#' the missing value in the indicated column.
#'
#' @param data A dataframe containing the column for which the values have to be imputed.
#' @param Column A character specifying the name of the column for which the missing values
#' have to be imputed.
#'
#' @return A dataframe that contains imputed values for the column that had missing values
#' initially, or the original dataframe if the indicated column did not contain missing
#' values.
#' @export
#'
#' @examples
#' dat_birds <- read.csv('./data-raw/Test_european_birds.csv')
#' sub_1 <- droplevels(subset(dat_birds, ID == 1))
#' test <- impute_ma(data = sub_1, column = 'Trait_mean')
#' sub_21 <- droplevels(subset(dat_birds, ID == 21))
#' test21 <- impute_ma(data = sub_21, column = 'Trait_mean')
impute_ma <- function(data, column){

  onlycol <- data[column]
  if (any(is.na(onlycol))) {
    forMean <- cbind(onlycol[[column]][which(is.na(onlycol)) - 1],
                     onlycol[[column]][which(is.na(onlycol)) + 1])

    data[[column]][is.na(onlycol)] <- rowMeans(forMean)
    return(data)
  } else {
    return(data)
  }
}


#' Impute a median of all the values in the column
#'
#' \code{imput_median} imputes the median of all the values in the indicated column.
#'
#' @param data A dataframe containing the column for which the values have to be imputed.
#' @param Column A character specifying the name of the column for which the missing values
#' have to be imputed.
#'
#' @return A dataframe that contains imputed values for the column that had missing values
#' initially, or the original dataframe if the indicated column did not contain missing
#' values.
#' @export
#'
#' @examples
#' dat_birds <- read.csv('./data-raw/Test_european_birds.csv')
#' sub_1 <- droplevels(subset(dat_birds, ID == 1))
#' test <- impute_median(data = sub_1, column = 'Trait_SE')
#' sub_21 <- droplevels(subset(dat_birds, ID == 21))
#' test21 <- impute_median(data = sub_21, column = 'Trait_SE')
impute_median <- function(data, column){
  onlycol <- data[column]

  if (any(is.na(onlycol))) {
    data[[column]][is.na(onlycol)] <- median(data[[column]], na.rm = T)
    return(data)
  } else {
    return(data)
  }
}

#' Extract path coefficients from the output of \code{\link{fit_all_meta}} function
#'
#' \code{extr_coefs} extracts the path coefficients for a specified relation.
#'
#' @param obj An object returned by the  \code{\link{fit_all_meta}} function
#' @param Type_EfS Character specifying for which relation to extract the path coefficients,
#' These relations reflect different pathways in the fitted SEM, for example 'Demog_rate_mean<-det_Clim',
#' 'Demog_rate_mean<-Pop_mean', 'Demog_rate_mean<-Trait_mean'. For more details see \code{\link{fit_meta}}
#'
#' @return A dataframe that contains 26 columns including the required path coefficients and
#' the associated meta-data.
#' @export
#'
#' @examples  ## come up with a simpler example???
#' Coefs_Aut <- readRDS(file = './output_forSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#' meta_Phen_Surv_neg <- fit_all_meta(data_MA = Coefs_Aut_neg,
#' Clim = 'Temperature',
#' Demog_rate = 'Survival',
#' Trait_categ = 'Phenological',
#' Covar = NULL,
#' COV = 'Pvalue + WeathQ',
#' sel = 'Temp_Neg_Phen_Surv_Cov',
#' folder_name = './output_overall_sig/',
#' colr = c('black'),
#' optimize = rep('uobyqa', 11))
#' test <- extr_coefs(data = meta_Phen_Surv_neg, Relation = 'Tot_DemRate<-det_Clim')
extr_coefs <- function(obj, Type_EfS){
  sub <- obj$data_meta[[1]]$data_EfS[data$data_meta[[1]]$names == Type_EfS][[1]]
  return(sub)
}
