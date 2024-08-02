#' Fit all meta-analytical models for a specified subselection of
#' trait and climate using selected phylogeny
#'
#' \code{fit_all_acountphylo} fits all mixed-effects meta-analytical models
#' (per each relation in the SEM), to extract global effect sizes across
#' the studies in a specified subset of data as defined by trait, climate
#' and using a given phylogenetic tree
#'
#' @param vtree_folder Character specifying the folder where the phylogenetic
#' tree is to be read from. The tree is expected to be in Newick format and be
#' named as "vertX.tre", where X is an integer (from 1 to 100).
#' @param ind An integer specifying what phylogenetic tree to be read from the
#' folder given by vtree_folder.
#' @inheritParams fit_all_meta
#'
#' @export
#'
#' @return Returns a tibble that includes the estimated across-study effect sizes,
#' their standard errors, their significance and the AIC for each fitted mixed-effects
#' model (fitted with or without phylogeny accounted for by using phylogenetic variance-covariance
#' matrix). This tibble also includes variances for random effects.
#' @examples
#' Coefs_Aut_sp <- readRDS(file = './output_all_simpleSEM/PathCoefs_allMods_Temp_Weights_DD_Autocor_FilteredDur_WITH_Traits.rds')
#' phenT <- fit_all_acountphylo(data_MA = Coefs_Aut_sp,
#'                              vtree = "./data/phylogenies",
#'                              ind = 2, Trait_categ = "Phenological",
#'                              Clim = "Temperature")

fit_all_acountphylo <- function(data_MA, vtree_folder, ind,
                                Trait_categ = 'Phenological',
                                Clim = 'Temperature'){
  vert_tree <- ape::read.tree(paste0(vtree_folder, "/vert", ind, ".tre"))
  # tre_ult <- ape::compute.brlen(vert_tree, power = 1)

  ## update the species names to correspond to the ones on phylogeny
  Coefs_Aut <- data_MA %>%
    dplyr::mutate(Sp_phylo = dplyr::case_when(
      Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
      Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
      Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
      Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
      Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
      Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
      TRUE ~ Species)) %>%
    dplyr::filter(! Sp_phylo %in% c('Chrysemys picta', 'Chelonia mydas'))


  Coefs_Aut$Species <- unlist(lapply(1:nrow(Coefs_Aut), FUN = function(x){
    binary <- strsplit(as.character(Coefs_Aut$Sp_phylo[x]), " ")
    Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")
  }))

  # need a copy with the sp names, to be used to account for phylogeny (apart from that
  # we also account for among-sp variation by including species as random intercept
  # in the model)
  Coefs_Aut$Sp_phylo <- Coefs_Aut$Species

  # prepare the var-covar matrix based on the requested Trait_categ
  Coefs_sub <- subset(Coefs_Aut, Relation == 'Trait_mean<-det_Clim' &
                             Trait_Categ == Trait_categ)
  tre_sub <- ape::drop.tip(vert_tree, which(!vert_tree$tip.label %in% Coefs_sub$Species))

  Mat_sub <- ape::vcv.phylo(tre_sub, corr = TRUE)

  # fitting meta-analysis for cliimate-dependent relations with both covariates
  meta_Cov <- fit_all_meta(data_MA = Coefs_Aut,
                                Demog_rate = NULL,
                                Trait_categ = Trait_categ,
                                Clim = Clim,
                                Cov_fact = 'WeathQ',
                                COV = 'Pvalue',
                                sel = 'DoesNotMatter',
                                folder_name = NULL,  # we do not want to print these detailed res-ts for these X trees
                                colr = c('black'),
                                DD = 'n_effectGR',
                                simpleSEM = TRUE,
                                A = Mat_sub,
                                all_Relations = c('Trait_mean<-det_Clim',
                                                  'Ind_GR<-det_Clim',
                                                  'Tot_GR<-det_Clim'))
  # fitting meta-analysis for all other relations
  meta_other <- fit_all_meta(data_MA = Coefs_Aut,
                             Demog_rate = NULL,
                             Trait_categ = Trait_categ,
                             Clim = Clim,
                             Cov_fact = NULL,
                             COV = NULL,
                             sel = 'DoesNotMatter',
                             A = Mat_sub,
                             folder_name = NULL,
                             colr = c('black'),
                             DD = 'n_effectGR',
                             simpleSEM = TRUE,
                             all_Relations = c('GR<-det_Clim', 'GR<-Pop_mean',
                                               'GR<-Trait_mean'))
  colnames(meta_Cov$meta_res[[1]])[1] <- 'Variable'
  ef_all <- rbind(meta_other$meta_res[[1]], meta_Cov$meta_res[[1]]) %>%
    dplyr::mutate(., Trait_Categ = Trait_categ,
                  Clim = Clim, phylo = ind)
  return(ef_all)
}

