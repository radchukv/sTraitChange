#' Analyze the results of climate window analysis
#'
#' \code{analyse_climwin} analyses the results obtained
#' with climate window analysis (function \code{\link{climwin_proc}})
#' for the specified study ID
#'
#' @param biol_data Data frame with trait data for a given
#' population and species.
#' @param ID Numeric giving a unique ID of the current dataset for
#' a given population and species.
#' @param randwin Logical (TRUE/FALSE) specifying whether \code{\link[climwin]{randwin}}
#' was run together with \code{\link[climwin]{slidingwin}}.
#' @param out_clim Character specifying the library on the path where
#' the results of climate window analysis were stored.
#' @param out_for_SEM Character specifying the library on the path where
#' to store the results that are to be used as input for SEM.
#' @param metric Character specifying 'AIC' or 'C'. This define whether
#' a value of PDAIC or Pc will be returned. For more information
#' see the same option in the function \code{\link[climwin]{pvalue}}.
#' @param test_winDur A Boolean specifying whether to check if the identified
#' window lays within the specified range and the delta AIC value is larger than
#' the specified threshold.
#' @param deltaThresh A numeric specifying the minimum delta AIC with the null
#' model to accept a specific climatic window model.
#' @param MinDur A numeric specifying the minimum allowed window duration.
#' @param MaxDur A numeric specifying the maximum allowed window duration.
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
#'                           out_clim = 'output_climwin',
#'                           out_for_SEM = 'output_forSEM'
#'                           randwin = TRUE, metric = 'AIC',
#'                           MinDur = 1, MaxDur = 12,
#'                           oneGrid = TRUE, explanYear = TRUE,
#'                           endWindow = 12, RefMon = NA)

analyse_climwin <- function(ID, biol_data,
                            out_clim = 'output_climwin_test',
                            randwin = FALSE, metric = 'AIC',
                            MinDur = 1, MaxDur = 40,
                            deltaThresh = -7, test_winDur = FALSE,
                            out_for_SEM = 'output_SEM_test',
                            oneGrid = TRUE, explanYear = TRUE,
                            endWindow, RefMon = NA) {

  subs <- droplevels(biol_data[biol_data$ID == ID, ])
  if (randwin) {
    dat <- readRDS(paste0('./', out_clim, '/', subs$ID[1], '_',
                          subs$Species[1], '_', subs$Location[1],
                          '_', subs$Trait[1], '_OneGrid_', oneGrid,
                          '_explYear_', explanYear, '_EndWindow_',
                          endWindow,'_RefMon_', RefMon,
                          '_Rand', '.RDS'))
    climwin_out <- dat$climwin_output[[1]]
    biol <- dat$biol_data[[1]]
    biol_data_noNA <- biol %>%
      filter(!is.na(Trait_mean) & !is.na(Trait_SE))
    randwin_out <- dat$randwin_output[[1]]
    data_climwin <- climwin_out$Dataset
    data_climwin$WindowDur <- data_climwin$WindowOpen - data_climwin$WindowClose
    climdata <- dat$clim_data[[1]]
  } else {
    dat <- readRDS(paste0('./', out_clim, '/', subs$ID[1], '_',
                          subs$Species[1], '_', subs$Location[1],
                          '_', subs$Trait[1], '_OneGrid_', oneGrid,
                          '_explYear_', explanYear, '_EndWindow_',
                          endWindow, '_RefMon_', RefMon,'.RDS'))

    climwin_out <- dat$climwin_output[[1]]
    biol <- dat$biol_data[[1]]
    data_climwin <- climwin_out$Dataset
    data_climwin$WindowDur <- data_climwin$WindowOpen - data_climwin$WindowClose
    climdata <- dat$clim_data[[1]]
  }

  pdf(paste0('./', out_clim, '/', subs$ID[1], '_',
             subs$Species[1], '_', subs$Location[1],
             '_', subs$Trait[1], '_OneGrid_', oneGrid,
             '_explYear_', explanYear, '_EndWindow_',
             endWindow, '_RefMon_', RefMon, '_climwin.pdf'))
  par(mfrow = c(2,2))
  print(climwin::plotdelta(climwin_out$Dataset))
  # print(climwin::plotwin(climwin_out$Dataset))  ## seems like a bug: if the window is exactly 1 time step, does not work
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
                                 sample.size = nrow(climwin_out$BestModelData))
    if(class(pval_rand) == 'character') {
      pval_rand <- as.numeric(substr(pval_rand, 2, 6))
    }
    # for now outputting climate var as calculated with the best model
    if(pval_rand < 0.05){
      message('Results of slidingwin() are unlikely an issue of overfitting.\n',
              'Randomization pvalue = ', round(pval_rand, 5), '\n')
    } else {
      message('Results of slidingwin() are likely an issue of overfitting.\n',
              'Randomization pvalue = ', round(pval_rand, 5), '\n')
    }
    if(test_winDur){
      dat_out <- check_winDur(climwin_out = climwin_out,
                              clim = climdata,
                              biol_data = biol_data_noNA,
                              MinDur = MinDur, MaxDur = MaxDur,
                              deltaThresh = deltaThresh)
    } else {
      dat_out <- cbind(biol_data_noNA, Clim =
                         climwin_out$BestModelData[, c('climate')],
                       WinDur = data_climwin$WindowDur[1],
                       deltaAIC = round(data_climwin$deltaAICc[1], 2))
    }
    if(! is.null(dat_out)){
      res <- tibble::tibble(ID = biol$ID[1],
                            Species = biol$Species[1],
                            Location = biol$Location[1],
                            Country = biol$Country[1],
                            Trait = biol$Trait[1],
                            pvalue = pval_rand,
                            orig_biol = list(biol),
                            data_res = list(dat_out))
      # save these data for SEM
      saveRDS(object = res,
              file = paste0('./', out_for_SEM, '/', biol$ID[1], '_',
                            biol$Species[1], '_', biol$Location[1],
                            '_', biol$Trait[1], '_OneGrid_', oneGrid,
                            '_explYear_', explanYear, '_EndWindow_',
                            endWindow, '_RefMon_', RefMon,
                            '_ForSEM',  '.RDS'))
      return(res)
      }
  } else {
    message('No ramdomizations were supplied and therefore \n',
            'no randomization test can be conducted \n',
            'No output dataset created')
  }
}
