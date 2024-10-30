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
#' @importFrom magrittr "%>%"
#'
#' @examples
#' # prepare data
#' dataPaths <- dataPaths %>%
#'                   dplyr::mutate(Species = dplyr::case_when(
#'                          Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
#'                          Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
#'                          Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
#'                          Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
#'                          Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
#'                          Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
#'                          TRUE ~ Species))
#'
#' dataPaths$Species <- unlist(lapply(1:nrow(dataPaths), FUN = function(x){
#'   binary <- strsplit(as.character(dataPaths$Species[x]), " ")
#'   Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")}))
#' dataPaths$Sp_phylo <- dataPaths$Species
#'
#' # fit the models for: Trait_mean<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'
#' meta_Phen_Cov <- fit_all_meta(data_MA = dataPaths,
#'                               Demog_rate = NULL,
#'                               Trait_categ = 'Phenological',
#'                               Clim = 'Temperature',
#'                               Cov_fact = 'WeathQ',
#'                               COV = 'Pvalue',
#'                               sel = 'Temp_Phen_Cov',
#'                               folder_name = NULL,
#'                               colr = c('black', 'red'),
#'                               DD = 'n_effectGR',
#'                               simpleSEM = TRUE,
#'                               A = phyloMat,
#'                               all_Relations = c('Trait_mean<-det_Clim',
#'                               'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'))
#' test <- extr_coefs(obj = meta_Phen_Cov, Type_EfS = 'Tot_GR<-det_Clim')
extr_coefs <- function(obj, Type_EfS){
  sub <- obj$data_meta[[1]]$data_EfS[obj$data_meta[[1]]$names == Type_EfS][[1]]
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
                                                         !grepl(':', data$Parameter)])] <- stats_out$QM
    data$pval[which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                         !grepl(':', data$Parameter)])] <- stats_out$QMp
    data$DF[which(data$Parameter %in% data$Parameter[grepl(variable, data$Parameter) &
                                                       !grepl(':', data$Parameter)])] <- stats_out$m
    return(data)
  }
}

