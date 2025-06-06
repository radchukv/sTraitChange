#' @title data_biol
#' @description subset of raw time series data containing the
#' metadata for each study (e.g. authors, location, species)
#' and for each year the population mean value for the trait
#' and the population size
#' @format A data frame with 328 rows and 33 variables:
#' \describe{
#'   \item{\code{ID}}{integer Study ID}
#'   \item{\code{Study_Authors}}{factor Authors of the study}
#'   \item{\code{Journal}}{factor Journal}
#'   \item{\code{Year_pub}}{factor Year of publication}
#'   \item{\code{Title}}{factor Title of the paper}
#'   \item{\code{Duration}}{integer Duration of the study
#'   in years}
#'   \item{\code{Species}}{factor Latin species name}
#'   \item{\code{Taxon}}{factor Taxon (can be bird, reptile,
#'   fish or mammal)}
#'   \item{\code{BirdType}}{factor Specifies whether the
#'   bird is a Seabird or not}
#'   \item{\code{Location}}{factor Location of the study}
#'   \item{\code{Longitude}}{double Longitude}
#'   \item{\code{Latitude}}{double Latitude}
#'   \item{\code{Weather_data_quality}}{factor Specifies
#'   whether climate was extracted from the location
#'    or from a nearby fastland (see Methods for details)}
#'   \item{\code{Country}}{factor Country}
#'   \item{\code{Year}}{integer Year}
#'   \item{\code{Trait}}{factor Specific type of the trait
#'   measured in the study, as reported in the publication}
#'   \item{\code{Trait_Categ_det}}{factor Type of trait, whereby
#'   specific traits reported in the individual studies were
#'   grouped if they reflected similar traits}
#'   \item{\code{Trait_Categ}}{factor Trait category, can
#'   be either Morphological or Phenological}
#'   \item{\code{Unit_trait}}{factor Units in which
#'   traits were measured}
#'   \item{\code{Trait_mean}}{double Annual population mean
#'    trait value}
#'   \item{\code{Trait_SE}}{double Standard error of the
#'   annual trait as reported in the study}
#'   \item{\code{Record_date}}{factor Period when the trait was
#'   recorded}
#'   \item{\code{Record_time}}{factor Specifies whether the
#'   trait was recorded during the breeding season}
#'   \item{\code{Demog_rate}}{factor Demographic rate measured}
#'   \item{\code{Demog_rate_Categ}}{factor Demographic rate category,
#'   one of: Reproduction, Recruitment, Survival}
#'   \item{\code{Demog_rate_unit}}{factor Units of demographic rates}
#'   \item{\code{Demog_rate_mean}}{double Annual population
#'   mean demographic rate value}
#'   \item{\code{Demog_rate_SE}}{double Standard error of the
#'   annual demographic rate value}
#'   \item{\code{Pop}}{factor Type of population size reported}
#'   \item{\code{Count}}{factor Specifies whether the
#'   population size was obtained by counts ('Y') or not ('N')}
#'   \item{\code{Pop_mean}}{double Annual population size}
#'   \item{\code{Pop_SE}}{double Standard error of the annual
#'   population size, if reported in the study}
#'   \item{\code{Closed}}{factor Specifies whether the extraction
#'   of the data from the publication was finished}
#'}
"data_biol"
