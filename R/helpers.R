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


#' Replace the stats in the table by the correct values from the omnibus test
#'
#' \code{replace_stats} replaces the values in the table by the correct output of
#' omnibus test.
#'
#' @param data A dataframe containing the output of the model summary
#' @param variable A character specifying the name of the variable for which the omnibus
#' test is ot be computed.
#'
#' @return A dataframe that contains corrected statistics values for the specified variable.
#' @export
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
#' stats <- stats::coef(summary(mod_genLength))
#' Inter_trait <- stats::anova(mod, btt = grep(':Trait', stats$Parameter))
#' DR <- stats::anova(mod, btt = which(stats$Parameter %in% stats$Parameter[grepl('Demog_rate', stats$Parameter) & !grepl(':', stats$Parameter)]))
#' stats_new <- replace_stats(data = stats, variable = ':Trait', stats_out = Inter_trait)
#' stats_new1 <- replace_stats(data = stats, variable = 'Demog_rate', stats_out = DR)
replace_stats <- function(data, variable, stats_out){
  if(length(grep(':', variable)) > 0){
    data$Chi2[grep(variable, data$Parameter)[1]] <- stats_out$QM
    data$DF[grep(variable, data$Parameter)[1]] <- stats_out$m
    data$pval[grep(variable, data$Parameter)[1]] <- stats_out$QMp
    if(length(grep(variable, data$Parameter)) > 1){
      for(i in seq(from = 2, to = length(grep(variable, data$Parameter)))){
        data$Chi2[grep(variable, data$Parameter)[i]] <- NA
        data$DF[grep(variable, data$Parameter)[i]] <- NA
        data$pval[grep(variable, data$Parameter)[i]] <- NA
      }
    }
    return(data)
  } else {
    data$Chi2[which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                         !grepl(':', data$Parameter)])[1]] <- stats_out$QM
    data$pval[which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                         !grepl(':', data$Parameter)])[1]] <- stats_out$QMp
    data$DF[which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                       !grepl(':', data$Parameter)])[1]] <- stats_out$m

    if(length(which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                       !grepl(':', data$Parameter)])) > 1){
      for(i in seq(from = 2,
                   to = length(which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                                        !grepl(':', data$Parameter)])))){
        data$Chi2[which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                             !grepl(':', data$Parameter)])[i]] <- NA
        data$DF[which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                           !grepl(':', data$Parameter)])[i]] <- NA
        data$pval[which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                             !grepl(':', data$Parameter)])[i]] <- NA
        return(data)
      }
    }
    return(data)
  }
}

