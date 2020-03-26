#' Prepare a dataset for SEM based on the output of window analyses
#'
#' \code{prep_SEM_input} produces a dataset to be used for SEM based on
#' the results of window analyses and supplied study IDs
#'
#' @param prep_subset_climwin The object returned by \code{\link{prep_subset}}
#' and that was used for window analyses.
#' @param out_for_SEM Character specifying the library on the path where
#' to store the results that are to be used as input for SEM.
#' @param oneGrid Logical (TRUE/FALSE). Whether climatic data for window analyses
#' were extracted from a single grid cell into which the study location falls or
#' if the mean of five cells (the focal one and four neighbours) was used.
#' @param explanYear Logical (TRUE/FALSE). Whether year was included as an
#' explanatory variable in the baseline formula for window analysis.
#' @param endWindow A numeric specifying the furtherst number of time intervals
#' (set by cinterval) back from the reference day. See also the option range in
#' the function \code{\link[climwin]{slidingwin}}.
#' @param RefMon A numeric specifying the month in case the absolute window was
#' for all species. If NA was specified then a species-specific window was used for
#' each dataset based on the available data on phenology and the period when morphology
#' was recorded.
#' @param selIDs A vector specifying the study IDs to be selected to include in the
#' dataset for SEM analyses.
#'
#' @export
#'
#' @return Returns a dataframe that for each study ID provides the data needed to conduct
#' SEM, i.e. climatic data extracted with window analyses, mean yearly trait values (with
#' SEs if available), mean yearly demographic rate values (with SEs if available) and
#' yearly population size estimates (with SEs if available). For each study ID metadata (e.g.
#' the study species, location, coordinates) are also included.
#'
#' @examples
#' temp_eu_SEM <- prep_SEM_input(prep_subset_climwin = eu_noSea,
#' oneGrid = FALSE,
#' explanYear = TRUE,
#' endWindow = 104,
#' RefMon = NA,
#' selIDs = 152)

prep_SEM_input <- function(prep_subset_climwin,
                           out_for_SEM,
                           oneGrid = FALSE,
                           explanYear = TRUE,
                           endWindow = 104,
                           RefMon = NA,
                           selIDs = unique(prep_subset_climwin$Sel[[1]]$ID)){
  data_all <- NULL
  for(j in selIDs){
     # cat('id (J) checked is', j, '\n')

    # save biol_NY as a dataset
    biol_NY <- prep_subset_climwin$biol_NY[[1]]
    # select one study
    subs <- droplevels(biol_NY[biol_NY$ID == j, ])

    ## check whether the file exists
    fileName <- paste0('./', out_for_SEM, '/', subs$ID[1], '_',
                       subs$Species[1], '_', subs$Location[1],
                       '_', subs$Trait[1],
                       '_OneGrid_', oneGrid,
                       '_explYear_', explanYear,
                       '_EndWindow_', endWindow,
                       '_RefMon_', RefMon, '_ForSEM',  '.RDS')

    if(file.exists(fileName)){
      test_fSEM <- readRDS(fileName)

      dat <- test_fSEM$data_res[[1]]
      dat$NYears <- nrow(dat)
      data_all <- rbind(data_all, dat)
      all_ID <- sort(unique(biol_NY$ID))

      all_ID <- all_ID[all_ID != j]
      for(i in all_ID){
        sub_NY <- droplevels(subset(biol_NY, ID == i))
        ## first check whether it is a study that has all the NAs (for SE we have some studies like that)

        if (sum(is.na(sub_NY$Trait_SE)) != nrow(sub_NY)){
        sub_NY_noNA <- sub_NY %>%
          dplyr::filter(!is.na(Trait_mean) & !is.na(Trait_SE)) %>%
          dplyr::mutate(., NYears = dplyr::n())
        }else{
          sub_NY_noNA <- sub_NY %>%
            dplyr::filter(!is.na(Trait_mean)) %>%
            dplyr::mutate(., NYears = dplyr::n())
        }

         # cat('ID (i) checked is', i, '\n')
        if(unique(dat$Study_Authors) == unique(sub_NY_noNA$Study_Authors) &
           unique(dat$Species) == unique(sub_NY_noNA$Species) &
           unique(dat$Location) == unique(sub_NY_noNA$Location) &
           unique(dat$Trait) == unique(sub_NY_noNA$Trait) &
           unique(dat$NYears) == unique(sub_NY_noNA$NYears)){
          data_new_sub <- cbind(sub_NY_noNA[, c(1:33)], dat[, c(34:46)])  ## check this one
          data_all <- rbind(data_all, data_new_sub)
        }
      }
    }
  }

  return(data_all)
}


