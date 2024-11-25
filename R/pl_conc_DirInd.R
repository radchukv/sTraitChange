#' Plot of the direct effect of climate on population growth rate vs the
#' trait-mediated effect of climate on population growth rate
#'
#'
#' \code{pl_conc_DirInd} Plots for each study the direct effect of climate
#' on population growth rate vs the trait-mediated effect of climate on population
#' growth rate and overlays the across-study estimates obtained with meta-analyses
#'
#' @inheritParams plot_concept
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#'
#' @return A list with 1. a plot of study-specific estimates of direct effects
#' of climate on population growth rate (CG) vs the estimates of trait-mediated
#' effects of climate on population growth rate (CZG), overlayed with the across-study
#' effect sizes derived from the meta-analyses as well as the estimated relation between
#' CG and CZG using the meta-analysis; 2. correlation estimates between CZG and CG;
#' 3. fitted mixed-effects model(s) predicting CG by CZG, if required per sign of the
#' climatic effect on traits.
#'
#' @examples
#' # prepare the data to fit the model to extract indirect path
#' dataPaths_sp <- dataPaths %>%
#'                   dplyr::mutate(Species = dplyr::case_when(
#'                          Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
#'                          Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
#'                          Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
#'                          Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
#'                          Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
#'                          Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
#'                          TRUE ~ Species))
#'
#' dataPaths_sp$Species <- unlist(lapply(1:nrow(dataPaths_sp), FUN = function(x){
#'   binary <- strsplit(as.character(dataPaths_sp$Species[x]), " ")
#'   Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")}))
#' dataPaths_sp$Sp_phylo <- dataPaths_sp$Species
#'
#' # fit the models for: Trait_mean<-det_Clim', 'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'
#' meta_Phen_Cov <- fit_all_meta(data_MA = dataPaths_sp,
#'                               Demog_rate = NULL,
#'                               Trait_categ = 'Phenological',
#'                               Clim = 'Temperature',
#'                               Cov_fact = 'WeathQ',
#'                               COV = 'Pvalue',
#'                               sel = 'Temp_Phen_Cov',
#'                               folder_name = NULL,
#'                               colr = c('black', 'red'),
#'                               DD = 'n_effectGR',
#'                               simpleSEM = TRUE,
#'                               A = phyloMat,
#'                               all_Relations = c('Trait_mean<-det_Clim',
#'                               'Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'))
#' two <- dplyr::bind_rows(lapply(list('Ind_GR<-det_Clim', 'Tot_GR<-det_Clim'),
#' function(x){extr_coefs(obj = meta_Phen_Cov, Type_EfS = x)}))
#' # fixing inconsistent variable naming
#' Coefs_subst <- dataPaths_sp %>%
#' dplyr::rename(.,  SError = Std.Error)
#' two$P.Value = rep(NA, nrow(two))
#'
#' Coefs_subst <- Coefs_subst %>%
#'     dplyr::select(., -c(names(Coefs_subst)[! names(Coefs_subst) %in% names(two)]))
#' two <- two %>%
#'     dplyr::select(., -c(names(two)[! names(two) %in% names(Coefs_subst)]))
#' Coef_all <- rbind(Coefs_subst, two) %>%
#' dplyr::filter(Trait_Categ == 'Phenological')
#' # now use this dataset with all coefs as input for ES_dat
#' allES_T <- Coef_all %>%
#'  dplyr::mutate(Climatic_var = 'Temperature')
#'  wide_temp_all <- allES_T %>%
#'   dplyr::select(ID, Estimate, SError, P.Value, Relation) %>%
#'   tidyr::pivot_wider(id_cols = ID, names_from = Relation,
#'   values_from = c(Estimate, SError, P.Value), names_sep = '/')
#'
#'   metaD_temp_all <- allES_T %>%
#'   dplyr::distinct(ID, .keep_all = TRUE) %>%
#'   dplyr::select(-c(Estimate, SError, Relation,
#'                   P.Value, Pvalue, Count, Nyears,
#'                   WinDur, Ref.day, Ref.month,
#'                   WindowClose, deltaAIC, Trait_ageClass,
#'                   WeathQ, GenLength_y_IUCN))
#'
#' wide_tempES_all <- (merge(wide_temp_all, metaD_temp_all, by = 'ID'))
#'
#' # and estimating the global, across-study estimates
#' globES_T <-  meta_Phen_Cov$meta_res[[1]] %>%
#' dplyr::filter(Levels_Covar == 'intrcpt') %>%
#' dplyr::mutate(Trait_Categ = 'Phenological',
#'   REL = dplyr::case_when(Relation == 'Trait_mean<-det_Clim' ~ 'CZ',
#'                         Relation == 'GR<-Trait_mean' ~ 'ZG',
#'                         Relation == 'GR<-det_Clim' ~ 'CG',
#'                         Relation == 'Ind_GR<-det_Clim' ~ 'CZG',
#'                         Relation == 'Tot_GR<-det_Clim' ~ 'TotalCG'))
#' # plot
#' PhenT_CZGvsCG <- pl_conc_DirInd(Trait_categ = 'Phenological',
#'                       GlobES_dat = globES_T,
#'                       ES_dat = wide_tempES_all,
#'                       xlab = 'Phenology-mediated effect <br> of temperature on G (CZ)',
#'                       ylab = 'Direct effect of temperature on G (CG)')$plot
pl_conc_DirInd <- function(Trait_categ = 'Phenological',
                           GlobES_dat,
                           ES_dat,
                           ylab = 'Direct effect of climate on GR (CG)',
                           xlab = 'Trait-mediated effect of climate on G (CZG)'){

  ES_dat <- ES_dat %>%
    dplyr::filter(.data$Trait_Categ == Trait_categ)
  GlobES_dat <- GlobES_dat[GlobES_dat$REL %in% c('CZG', 'CG') &
                             GlobES_dat$Trait_Categ == Trait_categ, ]
  GlobES_dat %<>%
    dplyr::mutate(ltype = dplyr::case_when(pval_Covar <= 0.1 ~ '1',
                                           TRUE ~ '2'))

    cort <- stats::cor.test(ES_dat$`Estimate/Ind_GR<-det_Clim`, ES_dat$`Estimate/GR<-det_Clim`)
    corel <- broom::tidy(cort)

    `SError/Ind_GR<-det_Clim` <- NULL
    mod_ML <- metafor::rma.mv(`Estimate/GR<-det_Clim` ~ `Estimate/Ind_GR<-det_Clim`,
                                  V = `SError/Ind_GR<-det_Clim`^2,
                                  random = list(~ 1|Species, ~1|ID, ~1|Location),
                                  data = ES_dat,
                                  method = 'ML')
    ## data frame to add the slope estimated with the mixed-effects model -
    ## at the moment not used to not overcrowd the plto
    rib_dat <- data.frame(x = seq(min(ES_dat$`Estimate/Ind_GR<-det_Clim`),
                                        max(ES_dat$`Estimate/Ind_GR<-det_Clim`),
                                        length.out = 10))

    rib_dat %<>%
      dplyr::mutate(ymax = mod_ML$ci.ub[2] * .data$x,
                    ymin = mod_ML$ci.lb[2] * .data$x,
                    `Estimate/GR<-det_Clim` = 0)

    `Estimate/Ind_GR<-det_Clim` <- `Estimate/GR<-det_Clim` <- REL <- Estimate <- ltype <- NULL
    Rel_CZG_CG <- data.frame(slope = mod_ML$beta["`Estimate/Ind_GR<-det_Clim`", 1],
                             intercept = 0,
                             ltype = ifelse(mod_ML$pval[2] <= 0.1, '1', '2'))

    if (requireNamespace("ggtext", quietly = TRUE)) {
    pl_CZGvsCG <- ggplot2::ggplot(ES_dat,
                                  ggplot2::aes(x = `Estimate/Ind_GR<-det_Clim`,
                                               y = `Estimate/GR<-det_Clim`)) +
      ggplot2::geom_point(col = 'black', alpha = 0.45, cex = 4) +
      ggplot2::geom_hline(data = subset(GlobES_dat, REL == 'CG'),
                          ggplot2::aes(yintercept = Estimate,
                                       lty = ltype), col = 'black') +
      ggplot2::geom_vline(data = subset(GlobES_dat, REL == 'CZG'),
                          ggplot2::aes(xintercept = Estimate,
                                       lty = ltype), col = 'black') +
      # ggplot2::geom_abline(data = Rel_CZG_CG,
      #                     ggplot2::aes(slope = slope,
      #                                    intercept = intercept,
      #                                    lty = ltype),
      #             col = 'black') +
      # ggplot2::geom_ribbon(data = rib_dat,
      #           ggplot2::aes(x = x, ymin= ymin, ymax = ymax, fill = SignClEffect),
      #             alpha = 0.3) +
      ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = 'none',  ## to not interfere with other legends in the composite plot
                         strip.background = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         strip.text = ggplot2::element_text(size  =12),
                         panel.grid.major = ggplot2::element_blank(),
                         axis.title = ggplot2::element_text(size = 20),
                         axis.text = ggplot2::element_text(size = 15),
                         axis.title.x = ggtext::element_markdown(),
                         axis.title.y = ggtext::element_markdown()) +
      ggplot2::scale_linetype_manual(values = c('1' = 1,
                                       '2' = 2),
                            labels = c('1' = 'p <= 0.1',
                                       '2' = 'p > 0.1'))
    } else {
      message("to be able to produce this plot, you first must run install.packages('ggtext')!")
    }
    if(requireNamespace('ggExtra', quietly = TRUE)){
    fin_pl <- ggExtra::ggMarginal(pl_CZGvsCG, type="density",
                                  fill = 'black', col = 'black')
    } else {
      message("to be able to produce this plot, you first must run install.packages('ggExtra')!")
    }
    return(list(plot = fin_pl, corTest = corel, GLMM = mod_ML))

}
