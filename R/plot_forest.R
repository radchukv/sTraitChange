#' Forest plot of effect sizes per each study and global effect sizes across studies
#'
#' \code{plot_forest} Plots a forest plot displaying the effect sizes
#' for each study and global effects across the studies
#'
#' @param data_ES Data frame containing, for each study, the effect size estimates and their
#' standard errors, returned as 'data_EfS' by the function \code{\link{fit_meta_phylo}}.
#' @param data_globES Data frame containing the global effect sizes across the studies,
#' returned by function \code{\link{fit_meta_phylo}} as either 'data_Covar' or 'data' (depending on
#' whether the explanatory categorical variable was included in the meta-analytical model).
#' @param xlab Character specifying the label for the x axis.
#' @param colr Vector specifying the colours to be used for the data points. The length of the
#' vector should correspond to the number of the levels in the categorical explanatory variable
#' included in the meta-analytical model, or should be one (if the single global effect size
#' across all studies is to be plotted).
#' @param pdf_basename A character specifying the name of the .pdf for
#' the produced plot.
#' @param mar A vector specifying the plot margins, analogously to \code{\link[graphics]{par}}.
#' @param labels_ES Boolean specifying whether to plot the labels for all the single studies.
#' Defaults to TRUE.
#' @inheritParams fit_meta_phylo
#'
#' @export
#'
#' @return Plots a forest plot with effect sizes (and SEs) for each study,
#' and global effect size(s) on the bottom. The data used for the plot are
#' returned (invisibly).
#'
plot_forest <- function(data_ES,
                        data_globES,
                        Cov_fact = NULL, COV = NULL,
                        xlab = 'Effect of temperature on trait',
                        colr = c('blue', 'green', 'black'),
                        pdf_basename = NULL,
                        mar = c(4, 10, 2, 2),
                        labels_ES = TRUE) {

  ## prepare data for the plot
  ymax <- NROW(data_ES) + NROW(data_globES)
  ymin <- ymax - nrow(data_ES) + 1
  data_ES$y <- ymax:ymin

  data_ES$lwr <- data_ES$Estimate - stats::qnorm(0.025)*data_ES$SError
  data_ES$upr <- data_ES$Estimate + stats::qnorm(0.025)*data_ES$SError

  if(! is.null(Cov_fact)){
    if(is.null(COV)){
      len_Covar <- length(unique(data_ES[, Cov_fact]))
    } else {
      len_Covar <- length(unique(data_globES[, 'Levels_Covar']))
    }

    for(j in 1:len_Covar){
      data_ES$colour[data_ES[, Cov_fact] == levels(data_ES[, Cov_fact])[j]] <- colr[j]
    }
  } else {
    data_ES$colour <- colr
  }
  data_ES$label <- paste(data_ES$Species, data_ES$ID, data_ES$Trait, sep = '_')  ## for now this simplified label - later see if more sophisiticated is needed

  ## prepare the limits of the plot
  ylim_plot <- c(1, ymax)
  xlim_max <-  ceiling(max(c(abs(data_ES$lwr), abs(data_ES$upr)), na.rm = TRUE))
  xlim_plot <- c(-xlim_max, xlim_max)


  ## start pdf if name of file defined
  if (!is.null(pdf_basename)) {
    grDevices::pdf(file = paste0(pdf_basename, '.pdf'))
  }

  graphics::par(mar = mar)
  ## drawing the empty plot
  graphics::plot(NULL,
                 xlim = xlim_plot,
                 ylim = ylim_plot,
                 ylab = '',
                 yaxt = 'n',
                 xlab = xlab,
                 cex.axis = 1.2,
                 cex.lab = 1.3)

  ## adding middle lines
  graphics::abline(v = 0, lty = 3)
  graphics::abline(h = 0, lty = 1)

  ## adding data
  for (i in 1:nrow(data_ES)) {
    graphics::points(data_ES[i, 'y'] ~ data_ES[i, 'Estimate'],
                     pch = 4,
                     lwd = 1.4,
                     col = data_ES[i, 'colour'])

    graphics::arrows(x0 = data_ES[i, 'lwr'], x1 = data_ES[i, 'upr'],
                     y0 = data_ES[i, 'y'],
                     y1 = data_ES[i, 'y'],
                     code = 3,
                     length = 0.01,
                     angle = 90,
                     lwd = 1.4,
                     col = data_ES[i, 'colour'])

    # graphics::segments(x0 = min(xlim_plot) - 2,
    #                    x1 = data_ES[i, 'lwr'],
    #                    y0 = data_ES[i, 'y'],
    #                    y1 = data_ES[i, 'y'],
    #                    lty = 3,
    #                    col = data_ES[i, 'colour'])

    if(labels_ES){
      graphics::mtext(data_ES[i, 'label'], side = 2, line = 0.5, at = data_ES[i, 'y'],
                      las = 2, cex = 0.7, col = data_ES[i, 'colour'])
    }

  }

  miny_glob <- min(data_ES$y) - 1

  ## prepare dataframe with global effect sizes
  data_globES_prep <- data_globES %>%
    dplyr::rename(.data, lwr = .data$EfS_Low, upr = .data$EfS_Upper) %>%
    dplyr::mutate(.data, y = c(1:.data$miny_glob))

  if(! is.null(Cov_fact)){
    data_globES_prep <- data_globES_prep %>%
      tidyr::separate(.data, .data$Levels_Covar, into = c('Category', 'Level'),
                      sep = Cov_fact, fill = 'right')

    # if(is.null(COV)){
    # for (j in 1:nrow(data_globES_prep)){
    #   data_globES_prep$colour[j] <- unique(data_ES$colour[data_ES[, Cov_fact] == data_globES_prep$Level[j]])
    # }
    # if(nrow(data_globES_prep) > 1){
    #   data_globES_prep$label <- data_globES_prep$Level
    # }
    # } else {
      for(j in 1:nrow(data_globES_prep)){
      data_globES_prep$colour[j] <- unique(data_ES$colour[data_ES[, Cov_fact] == data_globES_prep$Level[j]])
      # }
      #colr[length(colr)]
      #data_globES_prep$label <- c(data_globES_prep$Level[! is.na(data_globES_prep$Level)],
      #                            data_globES_prep$Category[! (data_globES_prep$Category == '')])

      data_globES_prep$label[j] <- ifelse(! (is.na(data_globES_prep$Category[j]) | data_globES_prep$Category[j] == ''), data_globES_prep$Category[j],
                                          ifelse(!is.na(data_globES_prep$Level[j]), data_globES_prep$Level[j], 'NA'))
      }
    data_globES_prep$colour[is.na(data_globES_prep$colour)] <- "black"
  } else {
    if(nrow(data_globES_prep) > 1){
      data_globES_prep$label <- data_globES_prep$Variable
      data_globES_prep$colour <- c(colr, rep('grey', nrow(data_globES_prep)-1))
    } else {
      data_globES_prep$label <- 'Global effect size'
    data_globES_prep$colour <- colr}
  }


  ## plot for the glob effect sizes

  for (i in 1:miny_glob){
    graphics::points(data_globES_prep[i, 'y'] ~ data_globES_prep[i, 'Estimate'],
                     pch = 4,
                     lwd = 1.4,
                     col = data_globES_prep[i, 'colour'],
                     cex = 1.2)

    graphics::arrows(x0 = data_globES_prep[i, 'lwr'], x1 = data_globES_prep[i, 'upr'],
                     y0 = data_globES_prep[i, 'y'],
                     y1 = data_globES_prep[i, 'y'],
                     code = 3,
                     length = 0.01,
                     angle = 90,
                     lwd = 1.8,
                     col = data_globES_prep[i, 'colour'])

    ## and here also labels
    graphics::mtext(data_globES_prep[i, 'label'], side = 2, line = 0.5, at = data_globES_prep[i, 'y'],
                    las = 2, cex = 0.9, col = data_globES_prep[i, 'colour'])

  }


  ## end pdf if name of file defined
  if (!is.null(pdf_basename)) {
    grDevices::dev.off()
    message(paste0('a pdf named', paste0(pdf_basename, '.pdf'), ' has been created and saved!'))
  }

  return(invisible(tibble::tibble(data_ES = list(data_ES),
                                  data_globES = list(data_globES_prep))))
}
