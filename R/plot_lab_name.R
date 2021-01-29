#' Providing labels for axes in forest plot and 'relation plots', as well as names for the pdfs
#'
#' \code{plot_lab_name} Prepares axes labels for forest plots and 'relation plots', as well as
#' names for the pdfs with plots
#'
#' @param Relation Character specifying for which type of the effect size a meta-analysis
#' was conducted. Possible are: 'Demog_rate_mean<-det_Clim', 'Demog_rate_mean<-Pop_mean',
#' 'Demog_rate_mean<-Trait_mean', 'GR<-Demog_rate_mean', 'GR<-det_Clim', 'GR<-Pop_mean',
#' 'Ind_DemRate<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_DemRate<-det_Clim',
#' 'Tot_GR<-det_Clim', and 'Trait_mean<-det_Clim'. For more details see \code{\link{fit_meta}}.
#' @param Demog_rate Character specifying the level of the demographic rate for which
#' analyses were conducted.
#' @param Trait_categ Character specifying the level of the trait, for which analyses were
#' conducted.
#' @param Clim Character specifying the level of the climatic variable, for which analyses
#' were conducted.
#' @inheritParams fit_meta
#'
#' @export
#'
#' @return A list with an x axis lable for the forest plot, x- and y-axes labels for 'relation plots'
#' and a nem for the pdf.
#'
#' @examples
#' lab_name <- plot_lab_name(Relation = 'GR<-det_Clim',
#'                           Cov_fact = NULL,
#'                           Demog_rate = 'survival',
#'                           Trait_categ = 'phenological',
#'                           Clim = 'temperature')
#'
plot_lab_name <- function(Relation = 'GR<-det_Clim',
                          Cov_fact = NULL,
                          Demog_rate = 'Survival',
                          Trait_categ = 'Phenological',
                          Clim = 'Temperature'){

  # to have suitable word for the axes labels
  if(Trait_categ == 'Phenological'){
    Trait_categ = 'Phenology'
  } else {
    if (Trait_categ == 'Morphological'){
      Trait_categ = 'Morphology'
    }
  }

  if(Relation == 'GR<-det_Clim'){
    xlab <- gen_specific_labs(vect = 'Direct effect of climate on pop GR', Clim = Clim, Trait_categ = Trait_categ, Demog_rate = Demog_rate)
    xlab_slopes <- gen_specific_labs(vect = 'Climate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Population GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)

  }
  if(Relation == 'GR<-Demog_rate_mean'){
    xlab <- gen_specific_labs(vect = 'Effect of demographic rate on pop GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Demographic rate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Population GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }
  if(Relation == 'GR<-Trait_mean'){
    xlab <- gen_specific_labs(vect = 'Effect of trait on pop GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Trait', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Population GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }
  if(Relation == 'GR<-Pop_mean'){
    xlab <- gen_specific_labs(vect = 'Effect of population size on pop GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Population size', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Population GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }
  if(Relation == 'Demog_rate_mean<-det_Clim'){
    xlab <- gen_specific_labs(vect = 'Direct effect of climate on demographic rate',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Climate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Demographic rate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }
  if(Relation == 'Demog_rate_mean<-Trait_mean'){
    xlab <- gen_specific_labs(vect =  'Effect of trait on demographic rate',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Trait', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Demographic rate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }
  if(Relation == 'Demog_rate_mean<-Pop_mean'){
    xlab <- gen_specific_labs(vect =  'Effect of population size on demographic rate',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Population size', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Demographic rate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }
  if(Relation == 'Trait_mean<-det_Clim'){
    xlab <- gen_specific_labs(vect = 'Effect of climate on trait',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Climate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Trait', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }
  if(Relation == 'Ind_DemRate<-det_Clim'){
    xlab <- gen_specific_labs(vect = 'Indirect effect of climate on demographic rate',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Climate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Demographic rate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }

  if(Relation == 'Ind_GR<-det_Clim'){
    xlab <- gen_specific_labs(vect = 'Indirect effect of climate on pop GR',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Climate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Population GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }

  if(Relation == 'Tot_DemRate<-det_Clim'){
    xlab <- gen_specific_labs(vect = 'Total effect of climate on demographic rate',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Climate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Demographic rate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }

  if(Relation == 'Tot_GR<-det_Clim'){
    xlab <- gen_specific_labs(vect = 'Total effect of climate on pop GR',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Climate', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Population GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }

  if(Relation == 'Ind_GR<-Pop_mean'){
    xlab <- gen_specific_labs(vect = 'Indirect effect of population size on pop GR',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Population size', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Population GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }

  if(Relation == 'Tot_GR<-Pop_mean'){
    xlab <- gen_specific_labs(vect = 'Total effect of population size on pop GR',
                              Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    xlab_slopes <- gen_specific_labs(vect = 'Population size', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
    ylab_slopes <- gen_specific_labs(vect = 'Population GR', Clim = Clim, Demog_rate = Demog_rate, Trait_categ = Trait_categ)
  }

  if(! is.null(Cov_fact)){
    pdf_name <- paste0(generate_pdf_name(Relation = Relation,
                                         Demog_rate = Demog_rate,
                                         Trait_categ = Trait_categ,
                                         Clim = Clim), '_by_', Cov_fact)
  } else {
    pdf_name <- generate_pdf_name(Relation = Relation, Demog_rate = Demog_rate, Trait_categ = Trait_categ, Clim = Clim)
  }
  return(list(xlab = xlab, pref_pdf = pdf_name,
              xlab_slope = xlab_slopes, ylab_slope = ylab_slopes))
}



#' Turn general axes labels into specific ones
#'
#' \code{gen_specific_labs} prepares the axes labels for the forest
#' and 'relation' plots. It replaces general names in axes (e.g. demographic rate)
#' by specific ones (e.g. survival). This is an internal function
#' used by \code{\link{plot_forest}}.
#' @param vect character specifying the general axis label that is to be converted
#' into a specific one.
#' @inheritParams plot_lab_name
#'
#' @return A vector of length one with the axis label.
#'
#' @export
#'
gen_specific_labs <- function(vect,
                              Clim = 'Temperature',
                              Demog_rate = 'Survival',
                              Trait_categ = 'Phenology'){
  if(length(unlist(gregexpr(" ", vect))) > 1) {  # to accommodate for 'demographi rate'
    converted <- sub(pattern = 'trait', replacement = tolower(Trait_categ),
                     x = sub(pattern = 'climate', replacement = tolower(Clim),
                             x = sub(pattern = 'demographic rate', replacement = tolower(Demog_rate),
                                     x = vect, ignore.case = TRUE),
                             ignore.case = TRUE), ignore.case = TRUE)
  } else {
    converted <- sub(pattern = 'trait', replacement = Trait_categ,
                      x = sub(pattern = 'climate', replacement = Clim,
                              x = sub(pattern = 'demographic rate', replacement = Demog_rate,
                                      x = vect, ignore.case = TRUE),
                              ignore.case = TRUE), ignore.case = TRUE)
  }
  return(converted)
}


#' Generate a pdf name for the plot
#'
#' \code{generate_pdf_name} prepares the pdf name for the forest
#' and 'relation' plots. This is an internal function
#' used by \code{\link{plot_forest}}.
#' @inheritParams plot_lab_name
#'
#' @return A vector of length one with the pdf name.
#'
#' @export
#'
generate_pdf_name <- function(Relation = 'GR<-det_Clim',
                              Clim = 'Temperature',
                              Demog_rate = 'Survival',
                              Trait_categ = 'Phenology'){

pdf_name <- gsub(pattern = "_|<|-", replacement = "",
  x = sub(pattern = 'Pop_mean', replacement = 'Pop',
            x = sub(pattern = 'Trait_mean',
                    replacement = substr(Trait_categ, 1, 4),
            x = gsub(pattern = "Demog_rate_mean|DemRate",
                    replacement = substr(Demog_rate, 1, 4),
            x = sub(pattern = 'det_Clim',
                    replacement = substr(Clim, 1, 4),
                    x = Relation, ignore.case = TRUE), ignore.case = TRUE),
            ignore.case = TRUE), ignore.case = TRUE))
return(pdf_name)
}
