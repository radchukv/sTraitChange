#' Plot of predictions using the full model from non-linearity test
#'
#' \code{plot_nonLin} Plots predictions with the fixed part of the model (in black)
#' and with the random model structure (in grey).
#'
#' @param data Data frame with the new data that will be used to predict with the model,
#' should have all the variables used as either fixed or random effects in the model.
#' @param mod The model with which the predictions are to be made.
#' @param pdf_basename A character specifying the name of the .pdf for
#' the produced plot.
#'
#' @export
#'
#' @return Plots the predictions and saves a pdf plot, if asked to.
#'
#' ## still add examples
#'
plot_nonLin <- function(data = NewDat_phen,
                        mod = Model,
                        pdf_basename = NULL){
  data$det_Clim2 <- data$det_Clim ^2
  data$Trait_mean2 <- data$Trait_mean ^2

  data <- data[order(data$ID, data$Year, data$det_Clim, data$Trait_mean), ]

  data$GR_pred <- predict(mod, newdata = data)
  data$GR_pred_fix <- predict(mod, newdata = data, re.form = NA) ## only fixed effects
  data$GR_pred_fix <- as.numeric(data$GR_pred_fix)
  data$GR_pred <- as.numeric(data$GR_pred)



  ## start pdf if name of file defined
  if (!is.null(pdf_basename)) {
    grDevices::pdf(file = paste0(pdf_basename, '.pdf'))
  }

  plot_all <- ggplot2::ggplot(data, aes(x = Trait_mean, y = GR_pred, group = ID)) +
    ggplot2::geom_point(colour = 'grey') +
    ggplot2::geom_line(data = data, aes(x = Trait_mean, y = GR_pred_fix),
                       col = 'black')  +
    ggplot2::geom_hline(yintercept = 0, col = 'darkgrey', lty = 3) +
    ggplot2::geom_vline(xintercept = 0, col = 'darkgrey', lty = 3) +
    xlab('Trait, Z | Climate, C') + ylab('Population growth rate, G') +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = 'bottom',
                                         strip.background = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         strip.text = element_text(size  =12),
                                         panel.grid.major = element_blank())
  print(plot_all)

  ## end pdf if name of file defined
  if (!is.null(pdf_basename)) {
    grDevices::dev.off()
    message(paste0('a pdf named', paste0(pdf_basename, '.pdf'), ' has been created and saved!'))
  }

}
