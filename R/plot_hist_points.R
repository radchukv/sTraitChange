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
#' @param Traitdem Character specifying the combination of the trait and demographic rate
#' categories, for which analyses were conducted. In case analyses were conducted separately
#' for trait category only, set to NULL.
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
#' @inheritParams fit_all_meta
#' @export
#'
#' @return Plots a forest-like plot with, in addition to the overall effect size estimates
#' shown by points and error bars, shown historgams of all the underlying study-speciifc
#' effect sizes. Each panel corresponds to a single relation (see \code{Type_EfS} in
#' \code{\link{fit_meta_phylo}} function) and effect sizes for precipitation and temperature
#' are shown in different colors and shapes.
#'
plot_hist_points <- function(data_allEstim,
                           Ef_sizesEstim,
                           dataAxis,
                           Traitdem = "Phenology",
                           tabC,
                           angle_yax_lab = 90, annot = NULL,
                           Trait_categ = NULL){
  Estimate <- Climatic_var <- EfS_Low <- EfS_Upper <-Col <- REL_Clim <- x <- y <- label <- Count <- NULL

  if(is.null(annot)){
    if(! is.null(Traitdem)){
      data_sub <- data_allEstim %>%
        dplyr::filter(.data$TraitDem == Traitdem)
      Ef_sub <- Ef_sizesEstim %>%
        dplyr::filter(.data$TraitDem == Traitdem)
      tab_subs <- tabC %>%
        dplyr::filter(.data$TraitDem == Traitdem)

    ggplot2::ggplot() + ggplot2::geom_histogram(data = data_sub,
                                                ggplot2::aes(x = Estimate, fill = Climatic_var),
                                                alpha = 0.5, bins = 50) +
        ggplot2::scale_fill_manual(values = c('Temperature' = 'lightsalmon1',
                                   'Precipitation' = 'lightblue1')) +
        ggplot2::geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', lwd = 0.9) +
        ggplot2::theme_bw() +
        ggplot2::theme(strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_blank(),
            strip.placement = "outside",
            panel.spacing = ggplot2::unit(0, "lines"),
            legend.position = 'bottom',
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(c(0.2, 0.2, 1.5, 1.5), unit = 'line'),
            legend.direction = 'vertical',
            panel.border = ggplot2::element_rect(fill = NA, color = 'grey',
                                        size = 0.1, linetype = 3))  +  ## , panel.border = element_blank()
        ggplot2::geom_errorbar(data = Ef_sub, width=.1,
                    ggplot2::aes(xmin = EfS_Low, xmax = EfS_Upper, colour = Col,
                                 y = 5), lwd = 0.5) +
        ggplot2::geom_point(data = Ef_sub,
                            ggplot2::aes(y = 5, x = Estimate, shape = Climatic_var,
                                         color = Col), size = 1.5) +
        ggplot2::facet_grid(rows = ggplot2::vars(REL_Clim)) +
        ggplot2::scale_colour_manual(values = c("Signif, temperature" = "red4",
                                     "Non-signif, temperature" = "tomato1",
                                     "Signif, precipitation" = "royalblue1",
                                     "Non-signif, precipitation" = "deepskyblue2")) +
        ggplot2::scale_shape_manual(values = c(16, 15)) +
        ggplot2::guides(color = ggplot2::guide_legend(ncol= 2, title = 'Overall effect size'),  ## ncol = 1 for ppt, where the legend is better on the right
             fill = ggplot2::guide_legend(ncol= 1, title = 'All effect sizes',
                                 override.aes = list(shape = c(NA, NA))),
             shape = ggplot2::guide_legend(title = 'Climatic variable')) +
        ggplot2::coord_cartesian(clip = 'off', xlim = c(-2, 2)) +
        ggplot2::geom_text(data = dataAxis, ggplot2::aes(x = x, y = y, label = label),
                           angle= angle_yax_lab, vjust = 0.5, hjust = 0.6) +
        ggplot2::labs(x = 'Effect size', y = '') + ggplot2::ggtitle(Traitdem) +
        ggplot2::geom_text(data = tab_subs,
                           ggplot2::aes(x = x, y = y, label = Count),
                           col = c('red4', 'royalblue1'))
    } else {
      data_sub <- data_allEstim %>%
        dplyr::filter(.data$Trait_Categ == Trait_categ)
      Ef_sub <- Ef_sizesEstim %>%
        dplyr::filter(.data$Trait_Categ == Trait_categ)
      tab_subs <- tabC %>%
        dplyr::filter(.data$Trait_Categ == Trait_categ)
      ggplot2::ggplot() + ggplot2::geom_histogram(data = data_sub,
                                                  ggplot2::aes(x = Estimate, fill = Climatic_var),
                                                  alpha = 0.5, bins = 50) +
        ggplot2::scale_fill_manual(values = c('Temperature' = 'lightsalmon1',
                                     'Precipitation' = 'lightblue1')) +
        ggplot2::geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', lwd = 0.9) +
        ggplot2::theme_bw() +
        ggplot2::theme(strip.background = ggplot2::element_blank(),
              strip.text = ggplot2::element_blank(),
              strip.placement = "outside",
              panel.spacing = ggplot2::unit(0, "lines"),
              legend.position = 'bottom',
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              plot.margin = ggplot2::margin(c(0.2, 0.2, 1.5, 1.5), unit = 'line'),
              legend.direction = 'vertical',
              panel.border = ggplot2::element_rect(fill = NA, color = 'grey',
                                          size = 0.1, linetype = 3))  +  ## , panel.border = element_blank()
        ggplot2::geom_errorbar(data =Ef_sub,
                      width=.1,
                      ggplot2::aes(xmin = EfS_Low, xmax = EfS_Upper, colour = Col, y = 5), lwd = 0.5) +
        ggplot2::geom_point(data = Ef_sub,
                   ggplot2::aes(y = 5, x = Estimate, shape = Climatic_var, color = Col), size = 1.5) +
        ggplot2::facet_grid(rows = ggplot2::vars(REL_Clim)) +
        ggplot2::scale_colour_manual(values = c("Signif, temperature" = "red4",
                                       "Non-signif, temperature" = "tomato1",
                                       "Signif, precipitation" = "royalblue1",
                                       "Non-signif, precipitation" = "deepskyblue2")) +
        ggplot2::scale_shape_manual(values = c(16, 15)) +
        ggplot2::guides(color = ggplot2::guide_legend(ncol= 2, title = 'Overall effect size'),  ## ncol = 1 for ppt, where the legend is better on the right
               fill = ggplot2::guide_legend(ncol= 1, title = 'All effect sizes',
                                   override.aes = list(shape = c(NA, NA))),
               shape = ggplot2::guide_legend(title = 'Climatic variable')) +
        ggplot2::coord_cartesian(clip = 'off', xlim = c(-2, 2)) +
        ggplot2::geom_text(data = dataAxis,
                           ggplot2::aes(x = x, y = y, label = label),
                           angle= angle_yax_lab, vjust = 0.5, hjust = 0.6) +
        ggplot2::labs(x = 'Effect size', y = '') + ggplot2::ggtitle(Trait_categ) +
        ggplot2::geom_text(data = tab_subs,
                           ggplot2::aes(x = x, y = y, label = Count),
                           col = c('red4', 'royalblue1'))
    }
  } else {
    hue <- NULL
    data_sub <- data_allEstim %>%
      dplyr::filter(.data$TraitDem == Traitdem)
    Ef_sub <- Ef_sizesEstim %>%
      dplyr::filter(.data$TraitDem == Traitdem)
    tab_subs <- tabC %>%
      dplyr::filter(.data$TraitDem == Traitdem)
 if(requireNamespace('ggnewscale', quietly = TRUE)){
    ggplot2::ggplot() + ggplot2::geom_blank(data =data_sub,
                                            ggplot2::aes(x = Estimate)) +
      ggplot2::geom_rect(data = annot,
                         ggplot2::aes(fill = hue), xmin = -Inf, xmax = Inf,
                         ymin = -Inf, ymax = Inf) +
      ggplot2::scale_fill_manual(values = c('Complex SEM' = 'lightgrey',
                                   'Simple SEM' = 'white'),
                        guide =  ggplot2::guide_legend(title = 'Model')) +
      ggnewscale::new_scale_fill() +
      ggplot2::geom_histogram(data = data_sub,
                              ggplot2::aes(x = Estimate, fill = Climatic_var), alpha = 0.5, bins = 50) +
      ggplot2::scale_fill_manual(values = c('Temperature' = 'lightsalmon1',
                                   'Precipitation' = 'lightblue1')) +
      ggplot2::geom_vline(xintercept = 0, linetype = 'dashed', color = 'black', lwd = 0.9) +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_blank(),
            strip.text = ggplot2::element_blank(),
            strip.placement = "outside",
            panel.spacing = ggplot2::unit(0, "lines"),
            legend.position = 'bottom',
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(c(0.2, 0.2, 1.5, 1.5), unit = 'line'),
            legend.direction = 'vertical',
            panel.border = ggplot2::element_rect(fill = NA, color = 'grey',
                                        size = 0.1, linetype = 3))  +  ## , panel.border = element_blank()
      ggplot2::geom_errorbar(data = Ef_sub,
                    width=.1,
                    ggplot2::aes(xmin = EfS_Low, xmax = EfS_Upper,
                                 colour = Col, y = 5), lwd = 0.5) +
      ggplot2::geom_point(data = Ef_sub, ggplot2::aes(y = 5, x = Estimate,
                                                      shape = Climatic_var, color = Col), size = 1.5) +
      ggplot2::facet_grid(rows = ggplot2::vars(REL_Clim)) +
      ggplot2::scale_colour_manual(values = c("Signif, temperature" = "red4",
                                     "Non-signif, temperature" = "tomato1",
                                     "Signif, precipitation" = "royalblue1",
                                     "Non-signif, precipitation" = "deepskyblue2")) +
      ggplot2::scale_shape_manual(values = c(16, 15)) +
      ggplot2::guides(color = ggplot2::guide_legend(ncol= 2,
                                                    title = 'Overall effect size'),  ## ncol = 1 for ppt, where the legend is better on the right
             fill = ggplot2::guide_legend(ncol= 1, title = 'All effect sizes',
                                 override.aes = list(shape = c(NA, NA))),
             shape = ggplot2::guide_legend(title = 'Climatic variable')) +
      ggplot2::coord_cartesian(clip = 'off', xlim = c(-2, 2)) +
      ggplot2::geom_text(data = dataAxis,
                         ggplot2::aes(x = x, y = y, label = label),
                         angle= angle_yax_lab,
                vjust = 0.5, hjust = 0.6) +
      ggplot2::labs(x = 'Effect size', y = '') + ggplot2::ggtitle(Traitdem) +
      ggplot2::geom_text(data = tab_subs,
                         ggplot2::aes(x = x, y = y, label = Count),
                         col = c('red4', 'royalblue1'))
 } else {
   message("to be able to produce this plot, you first must run install.packages('ggnewscale')!")
 }
  }
}
