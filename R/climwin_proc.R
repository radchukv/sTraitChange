#' Conduct climate window analysis
#'
#' \code{climwin_proc} performs climate window analysis for
#' specified traits and climate data and returns the best
#' identified window together with results of randomisation
#' tests (if asked for)
#'
#' @param biol_data Data frame with trait data for a given population and species.
#' @param clim_data Raster stack of climatic data across Europe
#' (if the study is outside Europe use XX).
#' @param ID Numeric giving a unique ID of the current dataset for a given population and species.
#' @param randwin Logical (TRUE/FALSE). Should \code{\link[climwin]{randwin}}
#'  be run together with \code{\link[climwin]{slidingwin}}?
#' @param repeats If randwin is run, the number of times that data will be randomised.
#' @param plot_check Logical (TRUE/FALSE). Whether to display
#' one year (the first in the series) of the climatic data
#'  together with the study location.
#' @param RefMon Numeric specifying the month for the refday(),
#' see the same option in the function \code{\link[climwin]{slidingwin}}
#' for more information.
#' @param out_dir Character specifying the library on the path where
#' to save the results
#' @param stat Character specifying which statistics to use for
#'  aggregating the climatic data (see also the same option in
#'  \code{\link[climwin]{slidingwin}}).
#'
#' @export
#'
#' @return A tibble with six columns if randwin is set to FALSE
#' and seven columns if randwin is set TRUE. The columns are: "ID"
#' - study id, "Species" - study species, "Trait" - the trait considered,
#' "climwin_output" - a list returned by slidingwin(),
#' "randwin_output" - a list returned by randwin (if set to TRUE),
#' "clim_data" - a data frame with climate values,
#' "biol_data" - a data frame with traits, demographic rates
#' and population size.
#'
#' @examples
#' dat_birds <- read.csv('./data-raw/Test_european_birds.csv')
#' meanT <- raster::stack('./data-raw/tg_ens_mean_0.1deg_reg_v18.0e.nc')
#' test_rand <- climwin_proc(biol_data = dat_birds,
#'                           clim_data = meanT, ID = 1,
#'                           randwin = FALSE, seednum = 1302,
#'                           repeats = 30, plot_check = FALSE,
#'                           RefMon = 6, out_dir = 'output_climwin',
#'                           stat = 'mean')
#'
climwin_proc <- function(biol_data, clim_data,
                         ID, randwin = FALSE,
                         seednum = 1302, repeats = 20,
                         plot_check = FALSE, RefMon = 6,
                         out_dir = 'output_climwin',
                         stat = 'mean'){

  biol_data <- droplevels(biol_data[biol_data$ID == ID, ])
  # add Date to biol data (for slidingwin)
  biol_data$Date <- as.Date(paste('01', '06', biol_data$Year,
                                  sep = '/'), format = '%d/%m/%Y')

  ## imputing mean of the previous and next values for the mean of the traits, and
  ## median of all the values for the SE of the traits
  biol_data <- impute_ma(data = biol_data, column = 'Trait_mean')
  biol_data <- impute_median(data = biol_data, column = 'Trait_SE')
  # biol_data <- biol_data[! is.na(biol_data$Trait_mean), ] -  former approach, excluding NAs

  location_spatial <-  sp::SpatialPointsDataFrame(coords = biol_data[1, c('Longitude', 'Latitude')],
                                                  data = data.frame(ID = biol_data$ID[1]),
                                                  proj4string = sp::CRS(as.character('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')))

  ##  for a visual check
  if (plot_check){
    raster::plot(clim_data[[1]])
    points(location_spatial)
  }

  ####                          Extract temperature data                  ####
  # Determine date data for each layer of the raster (allows us to sort by each year).
  temp_dates <- data.frame(Date = as.Date(names(clim_data),
                                          format = 'X%Y.%m.%d'))
  temp_dates$Year <- lubridate::year(temp_dates$Date)

  ## extract data from Euro weather for all the necessary dates for the site
  ## (i.e. everything from the year before the first recorded nest until the year of the most recent brood).
  message(paste('Currently extracting temperature data for study',
                biol_data$ID[1], 'for species',
                biol_data$Species[1], 'in', biol_data$Location[1], 'for',
                biol_data$Trait[1]))
  Clim <- data.frame(Date = seq(as.Date(paste('01', '01', min(biol_data$Year) - 1, sep = '/'), format = '%d/%m/%Y'),
                                as.Date(paste('01', '12', max(biol_data$Year), sep = '/'), format = '%d/%m/%Y'), 'day'),  # why max should be this and not + 1 year?
                     Temp = NA)  ##longer end date for clim dtaa compared to biol data is needed in order for basewin()
  ## checks to not crush

  ptstart <- proc.time()
  Clim$Temp <- ifelse(is.na(Clim$Temp),
                      as.numeric(raster::extract(clim_data[[which(temp_dates$Date %in%
                                                                    Clim$Date)]], location_spatial)),
                      NA)
  ptfinish <- proc.time() - ptstart
  cat('When extracting weather data from the raster the time elapsed is ',
      ptfinish[3], ';\n', 'time spend for ',
      names(ptfinish)[1], ' is ', ptfinish[1], ';\n',
      'for ', names(ptfinish)[2], ' is ', ptfinish[2], '\n',
      sep = '')

  # Check that data extracted properly!!
  if(all(is.na(Clim$Temp))){
    stop(paste0('Climate data failed to extract for ', biol_data$ID[1]))
  }

  # NOW THAT WE HAVE CLIMATE AND BIOLOGICAL DATA WE CAN RUN CLIMWIN!!
  set.seed(seednum)

  # create weights here, otherwise slidingwin does not see them
  biol_data$W <- 1 / biol_data$Trait_SE^2

  ptstart <- proc.time()
  climwin_output <- climwin::slidingwin(xvar = list(Temp = Clim$Temp),
                                        cdate = Clim$Date,
                                        bdate = biol_data$Date,
                                        baseline = lm(Trait_mean ~ 1, data = biol_data,
                                                      weights = W),
                                        range = c(365, 0),
                                        stat = stat, func = 'lin',
                                        type = 'absolute', refday = c(1, RefMon),
                                        cmissing = 'method2', cinterval = 'day')
  ptfinish <- proc.time() - ptstart
  cat('When fitting slidingwin() the time elapsed is ',
      ptfinish[3], ';\n', 'time spend for ',
      names(ptfinish)[1], ' is ', ptfinish[1], ';\n',
      'for ', names(ptfinish)[2], ' is ', ptfinish[2],
      '\n', sep = '')


  if(randwin){

    ptstart <- proc.time()
    randwin_output <- climwin::randwin(repeats = repeats,
                                       xvar = list(Temp = Clim$Temp),
                                       cdate = Clim$Date,
                                       bdate = biol_data$Date,
                                       baseline = lm(Trait_mean ~ 1, data = biol_data,
                                                     weights = W),
                                       range = c(365, 0),
                                       stat = stat, func = 'lin',
                                       type = 'absolute', refday = c(1, RefMon),
                                       cmissing = 'method2',
                                       cinterval = 'day')
    ptfinish <- proc.time() - ptstart
    cat('When running randwin() with ', repeats,
        'number of repeats, the time elapsed is ',
        ptfinish[3], ';\n', 'time spend for ',
        names(ptfinish)[1], ' is ', ptfinish[1], ';\n',
        'for ', names(ptfinish)[2], ' is ', ptfinish[2],
        '\n', sep = '')

    # create a tibble to save all output together
    clim_out <- tibble::tibble(ID = biol_data$ID[1],
                               Species = biol_data$Species[1],
                               Trait = biol_data$Trait[1],
                               climwin_output = list(climwin_output[[1]]),
                               randwin_output = list(randwin_output[[1]]),
                               clim_data = list(Clim),
                               biol_data = list(biol_data))
    saveRDS(object = clim_out,
            file = paste0('./', out_dir, '/', biol_data$ID[1], '_',
                          biol_data$Species[1], '_', biol_data$Location[1],
                          '_', biol_data$Trait[1], '_Rand',  '.RDS'))

  } else {
    # create a tibble to save all output together
    clim_out <- tibble::tibble(ID = biol_data$ID[1],
                               Species = biol_data$Species[1],
                               Trait = biol_data$Trait[1],
                               climwin_output = list(climwin_output[[1]]),
                               clim_data = list(Clim),
                               biol_data = list(biol_data))
    saveRDS(object = clim_out,
            file = paste0('./', out_dir, '/', biol_data$ID[1], '_',
                          biol_data$Species[1], '_',
                          biol_data$Location[1], '_', biol_data$Trait[1],
                          '.RDS'))
  }
  return(clim_out)
}
