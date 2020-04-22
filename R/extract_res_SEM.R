#' Extract results of the fitted SEM concerning a certain aspect of the model
#'
#' \code{extract_res_SEM} extracts the results of the list of models fitted with
#' \code{\link{fit_SEM}} function. The results are extracted with a focus on a particular aspect
#' of the model: Cstat value, R2 values or obtained path coefficients
#'
#' @param list_fitSEM List of SEMs fitted with the funciton \code{\link{fit_SEM}}.
#' @param stat_extr Character specifying on which aspect of the fitted SEM
#' you want to extract the results. Can take on the following values: 'coefs',
#' 'R2', 'Cstat'.
#'
#' @export
#'
#' @return A dataframe containing the meta-data on each study (e.g. study is,
#' species name, location, trait category, demographic rate category etc) and
#' the data on the requested aspect of the fitted SEM.
#'
#' still add the example

extract_res_SEM <- function(list_fitSEM, stat_extr){
  res <- dplyr::bind_rows(lapply(1:length(list_fitSEM), FUN = function(x){
    elem <- list_fitSEM[[x]]
    metadat <- elem[c('ID', 'Species', 'Location', 'Country',
                      'Continent', 'Longitude', 'Latitude',
                      'Taxon', 'BirdType', 'Trait_Categ', 'Trait',
                      'Demog_rate_Categ', 'Demog_rate', 'Count',
                      'Nyears', 'WinDur', 'deltaAIC', 'Pvalue',
                      'weights', 'DD', 'corr')]
    test <- elem[[stat_extr]][[1]]

    if(stat_extr == 'coefs'){
      test <- test[, -ncol(test)]
    }
    test <- test %>% dplyr::mutate_if(is.factor, as.character)
    ret_dat <- cbind(test,
                     metadat[rep(seq_len(nrow(metadat)), nrow(test)), ])
    # to avoid coersion to character when binding the rows by bind_rows
    ret_dat <- ret_dat %>% dplyr::mutate_if(is.factor, as.character)
  }
  ))
  return(res)
}
