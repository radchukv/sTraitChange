
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sTraitChange

<!-- badges: start -->
<!-- badges: end -->

The goal of sTraitChange is to provide the functions needed to reproduce
the analyses for the manuscript ‘Changes in phenology mediate vertebrate
population responses to climate globally’ by Radchuk et al. (submitted).

The package was conceived to perform analyses presented in the
above-mentioned paper. Therefore, the package will only be updated if
necessary to keep compatibility with other packages, so that the code
represented here works well.

## Installation

You can install the development version of sTraitChange from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("radchukv/sTraitChange")
```

## Usage

To access the documentation, type:

``` r
help(sTraitChange)
#> No documentation for 'sTraitChange' in specified packages and libraries:
#> you could try '??sTraitChange'
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
                   out_SEM = 'output_forSEM',
                   DD = 'n_effectGR', weight = TRUE,
                   correlation = TRUE,
                   standardize = TRUE,
                   Trait = FALSE,
                   simpleSEM = TRUE)
#> Currently fitting SEM for study 1 for species Parus major in Askainen for LayingDate
#> /private/var/folders/5z/w9py0h790j31_sv8r1xldf3w0000gq/T/Rtmpwe7mSw/temp_libpath119b05330e622/sTraitChange
#> /private/var/folders/5z/w9py0h790j31_sv8r1xldf3w0000gq/T/Rtmpwe7mSw/temp_libpath119b05330e622/sTraitChange/output_forSEM/1_Parus major_Askainen_LayingDate_z_score_relations.pdf
#> /private/var/folders/5z/w9py0h790j31_sv8r1xldf3w0000gq/T/Rtmpwe7mSw/temp_libpath119b05330e622/sTraitChange/output_forSEM/1_Parus major_Askainen_LayingDate_NumberFledglings_ResultsSEM.RDS
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

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

## Support

For any further information, please contact <radchuk@izw-berlin.de>.
