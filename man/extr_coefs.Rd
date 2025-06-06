% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{extr_coefs}
\alias{extr_coefs}
\title{Extract path coefficients from the output of \code{\link{fit_all_meta}} function}
\usage{
extr_coefs(obj, Type_EfS)
}
\arguments{
\item{obj}{An object returned by the  \code{\link{fit_all_meta}} function}

\item{Type_EfS}{Character specifying for which relation to extract the path coefficients,
These relations reflect different pathways in the fitted SEM, for example
'Demog_rate_mean<-det_Clim', Demog_rate_mean<-Pop_mean',
'Demog_rate_mean<-Trait_mean'. For more details see \code{\link{fit_meta_phylo}}}
}
\value{
A dataframe that contains 26 columns including the required path coefficients and
the associated meta-data.
}
\description{
\code{extr_coefs} extracts the path coefficients for a specified relation.
}
\examples{
# prepare data
dataPaths <- dataPaths \%>\%
                  dplyr::mutate(Species = dplyr::case_when(
                         Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
                         Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
                         Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
                         Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
                         Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
                         Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
                         TRUE ~ Species))

dataPaths$Species <- unlist(lapply(1:nrow(dataPaths), FUN = function(x){
  binary <- strsplit(as.character(dataPaths$Species[x]), " ")
  Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")}))
dataPaths$Sp_phylo <- dataPaths$Species

# fit the models for: Trait_mean<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'
meta_Phen_Cov <- fit_all_meta(data_MA = dataPaths,
                              Demog_rate = NULL,
                              Trait_categ = 'Phenological',
                              Clim = 'Temperature',
                              Cov_fact = 'WeathQ',
                              COV = 'Pvalue',
                              sel = 'Temp_Phen_Cov',
                              folder_name = NULL,
                              colr = c('black', 'red'),
                              DD = 'n_effectGR',
                              simpleSEM = TRUE,
                              A = phyloMat,
                              all_Relations = c('Trait_mean<-det_Clim',
                              'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'))
test <- extr_coefs(obj = meta_Phen_Cov, Type_EfS = 'Tot_GR<-det_Clim')
}
