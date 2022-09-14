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
#'
#' @return A list with 1. a plot of study-specific estimates of direct effects
#' of climate on population growth rate (CG) vs the estimates of trait-mediated
#' effects of climate on population growth rate (CZG), overlayed with the across-study
#' effect sizes derived from the meta-analyses as well as the estimated relation between
#' CG and CZG using the meta-analysis; 2. correlation estimates between CZG and CG;
#' 3. fitted mixed-effects model(s) predicting CG by CZG, if required per sign of the
#' climatic effect on traits.
#'
#' ## Still add examples
#'
pl_conc_DirInd <- function(Trait_categ = 'Phenological',
                           GlobES_dat = met_ef_T,
                           ES_dat = wide_tempES,
                           ClEfSpecific = TRUE){

  ES_dat <- subset(ES_dat, Trait_Categ == Trait_categ)
  GlobES_dat <- GlobES_dat[GlobES_dat$REL %in% c('CZG', 'CG') &
                             GlobES_dat$Trait_Categ == Trait_categ, ]
  GlobES_dat %<>%
    dplyr::mutate(ltype = dplyr::case_when(pval_across <= 0.1 ~ '1',
                                           TRUE ~ '2'))
  if(ClEfSpecific){
  corel <- ES_dat %>%
    tidyr::nest(data = -SignClEffect) %>%
    dplyr::mutate(
      cort = purrr::map(data, ~ cor.test(.x$`Estimate/Ind_GR<-det_Clim`,
                                         .x$`Estimate/GR<-det_Clim`)),
      out = purrr::map(cort, broom::tidy)) %>%
    tidyr::unnest(out)


  mod_ML_neg <- metafor::rma.mv(`Estimate/GR<-det_Clim` ~ `Estimate/Ind_GR<-det_Clim`,
                                V = `SError/Ind_GR<-det_Clim`^2,
                                random = list(~ 1|Species, ~1|ID, ~1|Location),
                                data = subset(ES_dat, SignClEffect == 'Negative'),
                                method = 'ML')
  mod_ML_pos <- metafor::rma.mv(`Estimate/GR<-det_Clim` ~ `Estimate/Ind_GR<-det_Clim`,
                                V = `SError/Ind_GR<-det_Clim`^2,
                                random = list(~ 1|Species, ~1|ID, ~1|Location),
                                data = subset(ES_dat, SignClEffect == 'Nonnegative'),
                                method = 'ML')

  ## data frame to add the slope estimated with the mixed-effects model -
  ## at the moment not used to not overcrowd the plto
  rib_dat <- data.frame(x = c(rep(seq(min(ES_dat$`Estimate/Ind_GR<-det_Clim`),
                                      max(ES_dat$`Estimate/Ind_GR<-det_Clim`),
                                      length.out = 10), 2)),
                        SignClEffect = rep(unique(ES_dat$SignClEffect), each = 10))

  rib_dat %<>%
    dplyr::mutate(ymax = c(mod_ML_neg$ci.ub[2] * x[SignClEffect == 'Negative'],
                           mod_ML_pos$ci.ub[2] * x[SignClEffect == 'Nonnegative']),
                  ymin = c(mod_ML_neg$ci.lb[2] * x[SignClEffect == 'Negative'],
                           mod_ML_pos$ci.lb[2] * x[SignClEffect == 'Nonnegative']),
                  `Estimate/GR<-det_Clim` = 0)

  Rel_CZG_CG <- data.frame(slope = c(mod_ML_neg$beta["`Estimate/Ind_GR<-det_Clim`", 1],
                                     mod_ML_pos$beta["`Estimate/Ind_GR<-det_Clim`", 1]),
                           SignClEffect = c('Negative', 'Nonnegative'), intercept = 0,
                           ltype = c(ifelse(mod_ML_neg$pval[2] <= 0.1, '1', '2'),
                                     ifelse(mod_ML_pos$pval[2] <= 0.1, '1', '2')))

  pl_CZGvsCG <- ggplot(ES_dat,
                       aes(x = `Estimate/Ind_GR<-det_Clim`,
                           y = `Estimate/GR<-det_Clim`,
                           col = SignClEffect)) +
    geom_point() +
    geom_hline(data = subset(GlobES_dat, REL == 'CG'),
               aes(yintercept = Estimate, col = SignClEffect,
                   lty = ltype)) +
    geom_vline(data = subset(GlobES_dat, REL == 'CZG'),
               aes(xintercept = Estimate, col = SignClEffect,
                   lty = ltype)) +
    geom_abline(data = Rel_CZG_CG, aes(slope = slope,
                                       intercept = intercept,
                                       col = SignClEffect,
                                       lty = ltype)) +
    # geom_ribbon(data = rib_dat, aes(x = x, ymin= ymin, ymax = ymax, fill = SignClEffect),
    #             alpha = 0.3) +
    xlab('Trait-mediated effect of climate on GR (CZG)') +
    ylab('Direct effect of climate on GR (CG)') +
    theme_bw() + theme(legend.position = 'bottom',
                       strip.background = element_blank(),
                       panel.grid.minor = element_blank(),
                       strip.text = element_text(size  =12),
                       panel.grid.major = element_blank()) +
    scale_color_manual(values = c('Negative' = 'darkorange',
                                  'Nonnegative' = 'darkgreen')) +
    scale_linetype_manual(values = c('1' = 1,
                                     '2' = 2),
                          labels = c('1' = 'p <= 0.1',
                                     '2' = 'p > 0.1')) +
    #geom_abline(intercept = 0, slope = -1, col = 'black') +
    guides(col = guide_legend(title = 'Climate effect'),
           lty = guide_legend(title = 'Significance'))

  fin_pl <- ggExtra::ggMarginal(pl_CZGvsCG, type="density",
                                groupFill = TRUE, groupColour = TRUE)
  return(list(plot = fin_pl, corTest = corel, GLMM_neg = mod_ML_neg,
              GLMM_pos = mod_ML_pos))
  } else {

    cort <- cor.test(ES_dat$`Estimate/Ind_GR<-det_Clim`, ES_dat$`Estimate/GR<-det_Clim`)
    corel <- broom::tidy(cort)

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
      dplyr::mutate(ymax = mod_ML$ci.ub[2] * x,
                    ymin = mod_ML$ci.lb[2] * x,
                    `Estimate/GR<-det_Clim` = 0)

    Rel_CZG_CG <- data.frame(slope = mod_ML$beta["`Estimate/Ind_GR<-det_Clim`", 1],
                             intercept = 0,
                             ltype = ifelse(mod_ML$pval[2] <= 0.1, '1', '2'))

    pl_CZGvsCG <- ggplot(ES_dat,
                         aes(x = `Estimate/Ind_GR<-det_Clim`,
                             y = `Estimate/GR<-det_Clim`)) +
      geom_point(col = 'black', alpha = 0.45) +
      geom_hline(data = subset(GlobES_dat, REL == 'CG'),
                 aes(yintercept = Estimate,  lty = ltype),
                 col = 'black') +
      geom_vline(data = subset(GlobES_dat, REL == 'CZG'),
                 aes(xintercept = Estimate, lty = ltype),
                 col = 'black') +
      geom_abline(data = Rel_CZG_CG, aes(slope = slope,
                                         intercept = intercept,
                                         lty = ltype),
                  col = 'black') +
      # geom_ribbon(data = rib_dat, aes(x = x, ymin= ymin, ymax = ymax, fill = SignClEffect),
      #             alpha = 0.3) +
      xlab('Trait-mediated effect of climate on GR (CZG)') +
      ylab('Direct effect of climate on GR (CG)') +
      theme_bw() + theme(legend.position = 'bottom',
                         strip.background = element_blank(),
                         panel.grid.minor = element_blank(),
                         strip.text = element_text(size  =12),
                         panel.grid.major = element_blank()) +
      scale_linetype_manual(values = c('1' = 1,
                                       '2' = 2),
                            labels = c('1' = 'p <= 0.1',
                                       '2' = 'p > 0.1')) +
      #geom_abline(intercept = 0, slope = -1, col = 'black') +
      guides(lty = guide_legend(title = 'Significance'))

    fin_pl <- ggExtra::ggMarginal(pl_CZGvsCG, type="density",
                                  fill = 'black', col = 'black')
    return(list(plot = fin_pl, corTest = corel, GLMM = mod_ML))
  }

}
