#' Plot of one chosen panel constituting conceptual plot (Z vs C, G vs Z|C or G vs C|Z)
#'
#' \code{plot_concept} Plots one chosen relation from the path diagram:
#' Z vs C, G vs Z|C or G vs C|Z for each specific study and adds the
#' across-study effect sizes estimated with meta-analytical models
#'
#' @param raw_dat A data frame with raw data that were used to fit single-study SEMs.
#' The data should be prepared in the way it was prepared for SEMs (e.g. climate, trait
#' and GR variables have to be standardized).
#' @param GlobES_dat A data frame containing the across-study effect size estimates
#' per relation, obtained with the meta-analyses fitted using the function
#' \code{\link{fit_all_meta}}.
#' @param ES_dat A data frame containing the study-specific effect size estimates per
#' relation (in a wide format), extracted from prop_data object obtained with the
#' function \code{\link{fit_all_meta}}.
#' @param path Character specifying which relation will be visualized on the plot.
#' Three possible: "CZ", "ZG" and "CG".
#' @param xvar_raw Character specifying what variable from the raw dataset (raw_dat)
#' plot on the x axis.
#' @param yvar_raw Character specifying what variable from the raw dataset (raw_dat)
#' plot on the y axis.
#' @param slope_ES Character specifying the name of the relation to be used to plot
#' from the ES_dat.
#' @param ylab Character specifying the label for the y axis.
#' @param ClEfSpecific Boolean indicating whether the meta-analyses were fitted separately for
#' studies with positive and negative effect of climate on traits.
#' @param miny Numeric specifying the minimum limit for the y axis.
#' @param maxy Numeric specifying the maximum limit for the y axis.
#'
#' @inheritParams fit_all_meta
#' @inheritParams plot_forest
#' @export
#'
#' @return Plots a requested relation (CZ, ZG or CG) for each study in the dataset and
#' overlays the across-study global effect sizes estimated with the meta-analyses.
#' ## Still add examples + revise to get the xvar_raw, yvar_raw, and slopeES programmatically!!
#'
plot_concept <- function(Trait_categ = 'Phenological',
                         raw_dat = temp_std,
                         GlobES_dat = met_ef_T,
                         ES_dat = wide_tempES,
                         path = 'CZ',
                         xvar_raw = 'det_Clim',
                         yvar_raw = 'Trait_mean',
                         slope_ES = 'Estimate/Trait_mean<-det_Clim',
                         ylab = 'Trait', xlab = 'Climate',
                         ClEfSpecific = TRUE,
                         miny = -6, maxy = 6){
  raw_dat <- subset(raw_dat, Trait_Categ == Trait_categ)
  GlobES_dat <- GlobES_dat[GlobES_dat$REL == path &
                             GlobES_dat$Trait_Categ == Trait_categ, ]
  GlobES_dat %<>%
    dplyr::mutate(ltype = dplyr::case_when(pval_Covar < 0.05 ~ '1',
                                           TRUE ~ '2'))
  ES_dat <- subset(ES_dat, Trait_Categ == Trait_categ)

  if(ClEfSpecific){
  dat_rib <- data.frame(x = c(rep(seq(min(raw_dat$det_Clim),
                                      max(raw_dat$det_Clim),
                                      length.out = 10), 2)),
                        SignClEffect = rep(unique(GlobES_dat$SignClEffect),
                                           each = 10))

  dat_rib %<>%
    dplyr::mutate(ymax = c(GlobES_dat$EfS_Upper[GlobES_dat$SignClEffect == 'Negative'] * x[SignClEffect == 'Negative'],
                           GlobES_dat$EfS_Upper[GlobES_dat$SignClEffect == 'Nonnegative'] * x[SignClEffect == 'Nonnegative']),
                  ymin = c(GlobES_dat$EfS_Low[GlobES_dat$SignClEffect == 'Negative'] * x[SignClEffect == 'Negative'],
                           GlobES_dat$EfS_Low[GlobES_dat$SignClEffect == 'Nonnegative'] * x[SignClEffect == 'Nonnegative']),
                  Trait_mean = 0, GR = 0)


  pl <- ggplot(raw_dat, aes(x = .data[[xvar_raw]],
                            y = .data[[yvar_raw]])) +
    lims(x = c(min(dat_rib$x), max(dat_rib$x)),
         y =  c(miny, maxy)) +
    geom_blank() +
    geom_abline(data = ES_dat,
                aes(intercept = 0, slope = .data[[slope_ES]]),
                col = 'lightgrey') +
    geom_abline(data = GlobES_dat, aes(intercept = 0, slope = Estimate,
                                       col = SignClEffect, lty = ltype),
                lwd = 1) +
    geom_ribbon(data = dat_rib, aes(x = x, ymin= ymin, ymax = ymax,
                                    fill = SignClEffect),
                alpha = 0.3) +
    scale_color_manual(values = c('Negative' = 'darkorange',
                                  'Nonnegative' = 'darkgreen')) +
    scale_linetype_manual(values = c('1' = 1,
                                     '2' = 2),
                          labels = c('1' = 'p <= 0.1',
                                     '2' = 'p > 0.1')) +
    theme_bw() + ylab(ylab) + xlab(xlab) +
    theme(legend.position = 'bottom',
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_markdown(),
          axis.title.y = element_markdown()) +
    guides(col = guide_legend(title = 'Climate effect'),
           fill = guide_legend(title = 'Climate effect'),
           lty = guide_legend(title = 'Significance'))
  } else {
    dat_rib <- data.frame(x = (seq(min(raw_dat$det_Clim),
                                        max(raw_dat$det_Clim),
                                        length.out = 10)))

    dat_rib %<>%
      dplyr::mutate(ymax = GlobES_dat$EfS_Upper * x,
                    ymin = GlobES_dat$EfS_Low * x,
                    Trait_mean = 0, GR = 0)


    pl <- ggplot(raw_dat, aes(x = .data[[xvar_raw]],
                              y = .data[[yvar_raw]])) +
      lims(x = c(min(dat_rib$x), max(dat_rib$x)),
           y =  c(miny, maxy)) +
      geom_blank() +
      geom_abline(data = ES_dat,
                  aes(intercept = 0, slope = .data[[slope_ES]]),
                  col = 'grey') +
      geom_ribbon(data = dat_rib, aes(x = x, ymin= ymin, ymax = ymax),
                  fill = 'black',
                  alpha = 0.55) +
      geom_abline(data = GlobES_dat,
                  aes(intercept = 0, slope = Estimate,
                      lty = ltype), col = 'black',
                  lwd = 1) +
      # scale_color_manual(values = c('Negative' = 'darkorange',
      #                               'Nonnegative' = 'darkgreen')) +
      scale_linetype_manual(values = c('1' = 1,
                                       '2' = 2),
                            labels = c('1' = 'p <= 0.1',
                                       '2' = 'p > 0.1')) +
      theme_bw() + ylab(ylab) + xlab(xlab) +
      theme(# legend.position = 'bottom',
            legend.position = 'none',  ## 22.06: to not interfere with other legends in the composite plot
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title = element_text(size = 25),
            axis.text = element_text(size = 20),
            axis.title.x = element_markdown(),
            axis.title.y = element_markdown()) +
      guides(lty = guide_legend(title = 'Significance'))
  }
  return(pl)
}
