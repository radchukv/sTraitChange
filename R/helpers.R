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
