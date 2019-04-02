#' Check whether the window duration is within realistic range
#'
#' \code{check_winDur} checks whether the window duration of the best
#' selected model is within the specified range and if not, finds the
#' model with such a window (given that the model with such window
#' is still better than the one without climatic variable as explanatory)
#'
#' @param biol_data Data frame with trait data for a given
#' population and species.
#' @param climwin_out A list containing the results of climatic window analyses,
#' returned by \code{\link[climwin]{slidingwin}} called within \code{\link{climwin_proc}}.
#' @param clim A dataframe consisting of daily climatic values,
#' corresponding to the clim_data element as returned by \code{\link{climwin_proc}}.
#' @param MinDur A numeric specifying the minimum allowed window duration.
#' @param MaxDur A numeric specifying the maximum allowed window duration.
#' @param deltaThresh A numeric specifying the minimum deltaAIC value for which the
#' model will stil be judged acceptable for selecting the window duration.
#'
#' @export
#'
#' @return If the climatic window within the specified range is found and the
#' respective model has deltaAIC at least -7 then the function returns a
#' dataframe consisting of biological data and an additional column 'Clim'
#' that for each year specifies the climatic values obtained according to
#' the selected model. If either the climatic window satisfying the specified range
#' is not found or / and the respective model has deltaAIC > -7, the function returns
#' NULL.
#'
#' @examples
#' dat_birds <- read.csv('./data-raw/Test_european_birds.csv')
#' meanT <- raster::stack('./data-raw/tg_ens_mean_0.1deg_reg_v18.0e.nc')
#' test_rand <- climwin_proc(biol_data = dat_birds,
#'                           clim_data = meanT, ID = 1,
#'                           randwin = TRUE, seednum = 1302,
#'                           repeats = 8, plot_check = FALSE,
#'                           RefMon = 7, out_dir = 'output_climwin_test',
#'                           stat = 'mean')
#' winDur <- check_winDur(climwin_out = test_rand$climwin_output[[1]],
#'                        clim = test_rand$clim_data[[1]],
#'                        biol_data = test_rand$biol_data[[1]],
#'                        MinDur = 7, MaxDur = 300, deltaThresh = -7)

check_winDur <- function(climwin_out, clim,
                         biol_data, MinDur,
                         MaxDur, deltaThresh = -7){
  data_climwin <- climwin_out$Dataset
  data_climwin$WindowDur <- data_climwin$WindowOpen - data_climwin$WindowClose

  if(data_climwin$WindowDur[1] > MinDur & data_climwin$WindowDur[1] < MaxDur){

    if(data_climwin$deltaAICc[1] <= deltaThresh){
      message('The best selected model with slidingwin() returns the', '\n',
              'climatic window within the specified range', '\n',
              'Window duration is ', data_climwin$WindowDur[1], '\n',
              'and the deltaAICc is satisfactory: ', round(data_climwin$deltaAICc[1], 2),  '\n')
      dat_out <- cbind(biol_data, Clim =
                         climwin_out$BestModelData[, c('climate')],
                       WinDur = data_climwin$WindowDur[1],
                       deltaAIC = round(data_climwin$deltaAICc[1], 2))
      return(dat_out)
    }else{
      message('The best selected model with slidingwin() returns the', '\n',
              'climatic window within the specified range', '\n',
              'Window duration is ', data_climwin$WindowDur[1], '\n',
              'but the deltaAICc is unsatisfactory: ', round(data_climwin$deltaAICc[1], 2),  '\n',
              'No output dataset created',  '\n')
    }

  }else{
    message('The best selected model with slidingwin() returns the', '\n',
            'climatic window outside the specified range', '\n',
            'Window duration is ', data_climwin$WindowDur[1], '\n')
    i <- 1
    while(data_climwin$WindowDur[i] < MinDur | data_climwin$WindowDur[i] > MaxDur){
      i <- i + 1
    }
    if(data_climwin$deltaAICc[i] <= deltaThresh){
      message('The identified climatic window within the specified range', '\n',
              'has the duration of ', data_climwin$WindowDur[i], '\n',
              'and the deltaAICc is satisfactory: ', round(data_climwin$deltaAICc[i], 2),  '\n')
      SelMod <- data_climwin[i, ]

      single_sel_win <- climwin::singlewin(xvar = list(Temp = clim$Temp),
                                           cdate = clim$Date,
                                           bdate = biol_data$Date,
                                           baseline = lm(Trait_mean ~ 1, data = biol_data,
                                                         weights = W),
                                           range = c(SelMod$WindowOpen, SelMod$WindowClose),
                                           stat = as.character(data_climwin$Statistics[1]),
                                           func = 'lin',
                                           type = 'absolute', refday = c(1, data_climwin$Reference.month[1]),
                                           cmissing = 'method2', cinterval = 'day')
      dat_out <- cbind(biol_data, Clim =
                         single_sel_win$BestModelData[, c('climate')],
                       WinDur = data_climwin$WindowDur[i],
                       deltaAIC = round(data_climwin$deltaAICc[i], 2))
      return(dat_out)
    } else {
      message('The identified climatic window within the specified range', '\n',
              'has the duration of ', data_climwin$WindowDur[i], '\n',
              'but the deltaAICc is unsatisfactory: ', round(data_climwin$deltaAICc[i], 2), '\n',
              'No output dataset created',  '\n')
    }
  }
}
