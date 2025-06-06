
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sTraitChange

<!-- badges: start -->
<!-- badges: end -->

The goal of sTraitChange is to provide the functions needed to reproduce
the analyses for the manuscript ‘Changes in phenology mediate vertebrate
population responses to temperature globally’ by Radchuk et
al. (submitted).

The package was conceived to perform analyses presented in the
above-mentioned paper. Therefore, the package will only be updated if
necessary to keep compatibility with other packages, so that the code
represented here works well.

## System requirements

### Hardware requirements

The `sTraitChange` package requires only a standard computer with enough
RAM to support the operations defined by a user. We recommend a computer
with following specs:  
- RAM: 16+ GB  
- CPU: 4+ cores.

### Software requirements

This package is supported for macOS and Windows operating systems. The
development version of `sTraitChange` was tested on:  
- macOS: Sequioa 15.0.1  
- Windows:

## Installation guide

The package depends on the following packages, that must be installed
prior to installing the `sTraitChange`:

``` r
install.packages(c('climwin', 'dplyr', 'tibble', 'ggplot2', 'magrittr', 'metafor', 'broom', 'piecewiseSEM', 'purrr', 'spaMM', 'tidyr', 'tidyselect', 'ape', 'sp', 'lubridate', 'openxlsx', 'ggnewscale', 'psych', 'ggtext', 'raster', 'ggExtra', 'data.table', 'rr2'))
```

The versions of these packages with which the package `sTraitChange` was
tested, are:

``` r
"ape 5.8"
"broom 1.0.7"
"climwin 1.2.3"
"data.table 1.16.2"
"dplyr 1.1.4"
"ggExtra 0.10.1"
"ggnewscale 0.5.0"
"ggplot2 3.5.1"
"ggtext 0.1.2"
"lubridate 1.9.3"
"magrittr 2.0.3"
"metafor 4.6.0"
"openxlsx 4.2.7.1"
"piecewiseSEM 2.3.0.1"
"psych 2.4.6.26"
"purrr 1.0.2"
"raster 3.6.30"
"sp 2.1.4"
"spaMM 4.5.0"
"tibble 3.2.1"
"tidyr 1.3.1"
"tidyselect 1.2.1"
```

You can install the development version of `sTraitChange` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("radchukv/sTraitChange")
```

This should take about 20 seconds on the computer with the recommended
specs.

## Demo

To access the documentation on a specific function, type ‘?’ and a
function name, for example:

``` r
?fit_SEM
```

The main functions of importance to the analyses performed for the paper
are:

- `climwin_proc()`: runs the sliding window analyses on each specific
  study;  
- `fit_SEM()`: fits the structural equation model to each study;  
- `fit_meta_phylo()`: fits the meta-analytical (mixed-effects model that
  accounts for phylogenetic relatedness) for one specific SEM path as a
  response variable;  
- `fit_all_meta()`: fits several meta-analytical models for several
  specified paths that will be used as response variables in the
  respective models.

And the main functions used for visualisation of the obtained results
are:  
- `plot_concept()`: plots the findings in an analogous way as
expectations laid out in the conceptual figure;  
- `plot_hist_points()`: plots the results of meta-analysis (i.e. global
effect sizes) overlaid over the histogram of the study-specific effect
sizes for specified paths from the SEM.

Here is an example of how SEM can be fitted to the 1st study in the
dataset:

``` r
mod_SEM <- fit_SEM(biol_data = dataSEM, ID = 1,
                   out_SEM = tempdir(), 
                   # attention: for this example we write the data to a temporary directory, to check its location type tempdir()
                   DD = 'n_effectGR', weight = TRUE,
                   correlation = TRUE,
                   standardize = TRUE,
                   Trait = FALSE,
                   simpleSEM = TRUE)
# check tempdir()
list.files(tempdir())
#> [1] "1_Parus major_Askainen_LayingDate_NumberFledglings_ResultsSEM.RDS"
#> [2] "1_Parus major_Askainen_LayingDate_relations.pdf"                  
#> [3] "1_Parus major_Askainen_LayingDate_z_score_relations.pdf"
```

And another example of how a meta-analysis that allows to account for
phylogenetic relatedness can be fitted to the phenological response to
temperature:

``` r
# prepare dataset, select only studies with phenological traits
Coefs_phenClim <- subset(dataPaths, 
                         Relation == 'Trait_mean<-det_Clim' & 
                           Trait_Categ =='Phenological')
# some standardization needed because phylogeny has different names for some species
Coefs_phenClim <- Coefs_phenClim %>%
                   dplyr::mutate(Species = dplyr::case_when(
                          Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
                          Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
                          Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
                          Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
                          Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
                          Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
                          TRUE ~ Species))
#  formatting sp names to have them exactly same way as on the phylogeny
Coefs_phenClim$Species <- unlist(lapply(1:nrow(Coefs_phenClim), FUN = function(x){
binary <- strsplit(as.character(Coefs_phenClim$Species[x]), " ")
Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")}))

# creating the object Sp_phylo that will account for the part in data variation
# purely due to phylogenetic relation
Coefs_phenClim$Sp_phylo <- Coefs_phenClim$Species

# fit meta-analysis
test_noCovar <- fit_meta_phylo(data_MA = Coefs_phenClim,
                               Type_EfS = 'Trait_mean<-det_Clim',
                               Cov_fact = NULL, COV = NULL,
                               DD = 'n_effectGR',
                               simpleSEM = TRUE,
                               A = phyloMat)
```

## Instructions for use

The applications of all the functions is demonstrated in detail in the
GitHub repo
[sTraitChange_Analyses](https://github.com/radchukv/sTraitChange_Analyses)
that implements the full workflow required to reproduce the results
presented in Radchuk et al. (submitted).  
Please check the help of the required functions and the repository
[sTraitChange_Analyses](https://github.com/radchukv/sTraitChange_Analyses)
for a specific type of the analyses you would like to apply on your
data.

## Support

For any further information, please contact <radchuk@izw-berlin.de>.
