#' Plot of the across-study effect sizes on top of the histogram
#' of all the effect sizes
#'
#' \code{plot_hist_points} Plots across-studies effect sizes as points and error bars on
#' top of the histogram of the underlying study-species effect sizes, per combination of
#' trait and demographic rate category
#'
#' @param data_allEstim A data frame with effect size estimates per relation for each
#' study (as returned by the SEM fitted with the function \code{\link{fit_SEM}}.
#' @param Ef_sizesEstim A data frame containing the across-study effect size estimates
#' per relation, obtained with the meta-analyses fitted using the function
#' \code{\link{fit_all_meta}}.
#' @param dataAxis A data frame specifying the names of the relations to be plotted on the
#' y axis.
#' @param Traitdem Character specifying the level of the trait, for which analyses were
#' conducted.
#' @param tabC A data frame containing the number of studies per each combination of
#' trait and demographic rate category.
#' @param angle_yax_lab Numeric specifying whether the labels of relations provided with
#' \code{dataAxis} are to be plotted parallel (\code{dataAxis} = 0) or perpendicular
#' (\code{dataAxis} = 90) to the x axis.
#' @param annot A data frame supplying the information on whether each relation was
#' estimated with complex or simple SEM (see \code{\link{fit_mod}} for more details),
#' so that the respective plot panesl will be shaded accordingly. Applicable only if
#' the plot combines the relations from simple and complex SEMs, otherwise set to NULL
#' (default).
#'
#' @export
#'
#' @return Plots a forest-like plot with, in addition to the overall effect size estimates
#' shown by points and error bars, shown historgams of all the underlying study-speciifc
#' effect sizes. Each panel corresponds to a single relation (see \code{Type_EfS} in
#' \code{\link{fit_meta}} function) and effect sizes for precipitation and temperature
#' are shown in different colors and shapes.
#' ## Still add examples
#'
plot_hist_points <- function(data_allEstim,
                           Ef_sizesEstim,
                           dataAxis,
                           Traitdem = "Morphology + Recruitment",
                           tabC,
                           angle_yax_lab = 90, annot = NULL){
  if(is.null(annot)){
    ggplot() + geom_histogram(data = subset(data_allEstim, TraitDem == Traitdem),
                              aes(x = Estimate, fill = Climatic_var), alpha = 0.5, bins = 50) +
      scale_fill_manual(values = c('Temperature' = 'lightsalmon1',
                                   'Precipitation' = 'lightblue1')) +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', lwd = 0.9) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            strip.placement = "outside",
            panel.spacing = unit(0, "lines"),
            legend.position = 'bottom',
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(c(0.2, 0.2, 1.5, 1.5), unit = 'line'),
            legend.direction = 'vertical',
            panel.border = element_rect(fill = NA, color = 'grey',
                                        size = 0.1, linetype = 3))  +  ## , panel.border = element_blank()
      geom_errorbar(data = subset(Ef_sizesEstim, TraitDem == Traitdem),
                    width=.1,
                    aes(xmin = EfS_Low, xmax = EfS_Upper, colour = Col, y = 5), lwd = 0.5) +
      geom_point(data = subset(Ef_sizesEstim, TraitDem == Traitdem),
                 aes(y = 5, x = Estimate, shape = Climatic_var, color = Col), size = 1.5) +
      facet_grid(rows = vars(REL_Clim)) +
      scale_colour_manual(values = c("Signif, temperature" = "red4",
                                     "Non-signif, temperature" = "tomato1",
                                     "Signif, precipitation" = "royalblue1",
                                     "Non-signif, precipitation" = "deepskyblue2")) +
      scale_shape_manual(values = c(16, 15)) +
      guides(color = guide_legend(ncol= 2, title = 'Overall effect size'),  ## ncol = 1 for ppt, where the legend is better on the right
             fill = guide_legend(ncol= 1, title = 'All effect sizes',
                                 override.aes = list(shape = c(NA, NA))),
             shape = guide_legend(title = 'Climatic variable')) +
      coord_cartesian(clip = 'off', xlim = c(-2, 2)) +
      geom_text(data = dataAxis, aes(x = x, y = y, label = label), angle= angle_yax_lab,
                vjust = 0.5, hjust = 0.6) +
      labs(x = 'Effect size', y = '') + ggtitle(Traitdem) +
      geom_text(data = subset(tabC, TraitDem == Traitdem),
                aes(x = x, y = y, label = Count), col = c('red4', 'royalblue1'))
  } else {
    ggplot() + geom_blank(data = subset(data_allEstim, TraitDem == Traitdem),
                          aes(x = Estimate)) +
      geom_rect(data = annot, aes(fill = hue),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      scale_fill_manual(values = c('Complex SEM' = 'lightgrey',
                                   'Simple SEM' = 'white'),
                        guide =  guide_legend(title = 'Model')) +
      ggnewscale::new_scale_fill() +
      geom_histogram(data = subset(data_allEstim, TraitDem == Traitdem),
                     aes(x = Estimate, fill = Climatic_var), alpha = 0.5, bins = 50) +
      scale_fill_manual(values = c('Temperature' = 'lightsalmon1',
                                   'Precipitation' = 'lightblue1')) +
      geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', lwd = 0.9) +
      theme_bw() +
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            strip.placement = "outside",
            panel.spacing = unit(0, "lines"),
            legend.position = 'bottom',
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(c(0.2, 0.2, 1.5, 1.5), unit = 'line'),
            legend.direction = 'vertical',
            panel.border = element_rect(fill = NA, color = 'grey',
                                        size = 0.1, linetype = 3))  +  ## , panel.border = element_blank()
      geom_errorbar(data = subset(Ef_sizesEstim, TraitDem == Traitdem),
                    width=.1,
                    aes(xmin = EfS_Low, xmax = EfS_Upper, colour = Col, y = 5), lwd = 0.5) +
      geom_point(data = subset(Ef_sizesEstim, TraitDem == Traitdem),
                 aes(y = 5, x = Estimate, shape = Climatic_var, color = Col), size = 1.5) +
      facet_grid(rows = vars(REL_Clim)) +
      scale_colour_manual(values = c("Signif, temperature" = "red4",
                                     "Non-signif, temperature" = "tomato1",
                                     "Signif, precipitation" = "royalblue1",
                                     "Non-signif, precipitation" = "deepskyblue2")) +
      scale_shape_manual(values = c(16, 15)) +
      guides(color = guide_legend(ncol= 2, title = 'Overall effect size'),  ## ncol = 1 for ppt, where the legend is better on the right
             fill = guide_legend(ncol= 1, title = 'All effect sizes',
                                 override.aes = list(shape = c(NA, NA))),
             shape = guide_legend(title = 'Climatic variable')) +
      coord_cartesian(clip = 'off', xlim = c(-2, 2)) +
      geom_text(data = dataAxis, aes(x = x, y = y, label = label), angle= angle_yax_lab,
                vjust = 0.5, hjust = 0.6) +
      labs(x = 'Effect size', y = '') + ggtitle(Traitdem) +
      geom_text(data = subset(tabC, TraitDem == Traitdem),
                aes(x = x, y = y, label = Count), col = c('red4', 'royalblue1'))
  }
}
