#' Analyze the results of climate window analysis
#'
#' \code{analyse_climwin} analyses the results obtained
#' with climate window analysis (function \code{\link{climwin_proc}})
#' for the specified study ID
#'
#' @param biol_data Data frame with trait data for a given population and species.
#' @param ID Numeric giving a unique ID of the current dataset for a given population and species.
#' @param randwin Logical (TRUE/FALSE) specifying whether \code{\link[climwin]{randwin}}
#' was run together with \code{\link[climwin]{slidingwin}}.
#' @param out_dir Character specifying the library on the path where
#' the results of climate window analysis were stored.
#' @param metric Character specifying 'AIC' or 'C'. This define whether
#' a value of PDAIC or Pc will be returned. For more information
#' see the same option in the function \code{\link[climwin]{pvalue}}
#'
#' @export
#'
#' @return If the results are not likely to be an issue of overfitting,
#' then the function returns a tibble with four columns: "ID"
#' - study id, "Species" - study species, "pvalue" - the pvalue from
#' the test assessing whether the climate signal is obtained by chance only,
#' and "data_res" - a data frame with all the information required
#' for running the SEM. If the results are likely an issue of
#' overfitting, then the returned tibble contains only first three of
#' the above-mentioned columns.
#'
#' @examples
#' dat_birds <- read.csv('./data-raw/Test_european_birds.csv')
#' t_anal <- analyse_climwin(ID = 1, biol_data = dat_birds,
#'                           out_dir = 'output_climwin',
#'                           randwin = TRUE, metric = 'C')

analyse_climwin <- function(ID, biol_data,
                            out_dir = 'output_climwin',
                            randwin = FALSE, metric = 'C'){

  subs <- droplevels(biol_data[biol_data$ID == ID, ])
  if (randwin) {
    dat <- readRDS(paste0('./', out_dir, '/', subs$Species[1], '_',
                          subs$Location[1], '_Rand',  '.RDS'))
    climwin_out <- dat$climwin_output[[1]]
    biol <- dat$biol_data[[1]]
    randwin_out <- dat$randwin_output[[1]]
  } else {
    dat <- readRDS(paste0('./', out_dir, '/', subs$Species[1], '_',
                          subs$Location[1], '.RDS'))

    climwin_out <- dat$climwin_output[[1]]
    biol <- dat$biol_data[[1]]
  }

  # climwin_out$BestModel
  # head(climwin_out$BestModelData)

  # discuss with Liam whether to take the results of the best
  # model or of the medwin
  # climwin::medwin(climwin_out$Dataset, cw = 0.95)

  pdf(paste0('./', out_dir, '/', subs$Species[1], '_',
             subs$Location[1], '_climwin.pdf'))
  par(mfrow = c(2,2))
  print(climwin::plotdelta(climwin_out$Dataset))
  print(climwin::plotwin(climwin_out$Dataset))
  print(climwin::plotweights(climwin_out$Dataset))
  if(randwin){
    print(climwin::plothist(dataset = climwin_out$Dataset,
                            datasetrand = randwin_out))
  }
  dev.off()

  if(randwin){
    pval_rand <- climwin::pvalue(dataset = climwin_out$Dataset,
                                 datasetrand = randwin_out,
                                 metric = metric,
                                 sample.size = nrow(biol))
    # for now outputting climate var as calculated with the best model
    if(pval_rand < 0.05){
      message('Results of slidingwin() are unlikely an issue of overfitting.\n',
              'Randomization pvalue = ', round(pval_rand, 5), '\n')
      dat_out <- cbind(biol, Clim =
                         climwin_out$BestModelData[, c('climate')])

      res <- tibble::tibble(ID = biol$ID[1],
                            Species = biol$Species[1],
                            pvalue = pval_rand,
                            data_res = list(dat_out))
      # save these data for SEM
      saveRDS(object = res,
              file = paste0('./', out_dir, '/', biol$Species[1], '_',
                            biol$Location[1], '_ForSEM',  '.RDS'))

      return(res)
    } else {
      message('Results of slidingwin() are likely an issue of overfitting.\n',
              'Randomization pvalue = ', round(pval_rand, 5), '\n',
              'No output dataset created')
      return(tibble::tibble(ID = biol$ID[1],
                            Species = biol$Species[1],
                            pvalue = pval_rand))
    }
  } else {
    message('No ramdomizations were supplied and therefore \n',
            'no randomization test can be conducted \n',
            'No output dataset created')
  }
}
