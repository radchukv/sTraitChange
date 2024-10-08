#' Plot predictions of climate effects on traits, using species-specific predictors
#'
#' \code{plot_uni_spSpec} Plots predictions of the model using as a response effect
#' of climate on trait and as predictor one species-specific predictor at a time and
#' climate
#'
#' @param data_allEstim A data frame with effect size estimates per relation for each
#' study (as returned by the SEM fitted with the function \code{\link{fit_SEM}}.
#' @param mod_mv Model object of class rma.mv fitted to test the
#' hypothesis about the effect of the specified species characteristic..
#' @param lOut Numeric specifying the length of the predictor variable, for which the
#' response of species traits to climate will be predicted and plotted.
#' @param pdf_basename A character specifying the name of the .pdf for
#' the produced plot.
#' @param xLab Character specifying the x axis label.
#' @param xLab Character specifying the y axis label.
#' @param byHemisphere Boolean specifying whether raw data points should be coloured depending on
#' which Hemisphere they are coming from.
#' @miny Numeric specifying the minimum limit for the y axis.
#' @maxy Numeric specifying the maximum limit for the y axis.
#' @export
#'
#' @return Plots a forest plot with effect sizes (and SEs) for each study,
#' and global effect size(s) on the bottom. The data used for the plot are
#' returned (invisibly).
#'
#' @examples
#' STILL To add
#'
plot_uni_spSpec <- function(data_allEstim = CZ_Phen,
                            mod_mv = mod_CZ_Lat,
                            lOut = 10,
                            pdf_basename = NULL,
                            xLab,
                            yLab = 'Estimate',
                            byHemisphere = FALSE,
                            miny = -1.2, maxy = 1.2){

  subs_pred <- names(as.data.frame(predict(mod_mv,  addx=TRUE)))[seq(length(names(as.data.frame(predict(mod_mv,  addx=TRUE)))) - nrow(stats::coef(summary(mod_mv))) + 1, length(names(as.data.frame(predict(mod_mv,  addx=TRUE)))))]
  subs_pred <- unlist(lapply(strsplit(subs_pred, split = 'X.'), function(x){x[2]}))
  Pred_data <- matrix(c(seq(from = min(data_allEstim[, subs_pred[2]], na.rm = T),
                                to = max(data_allEstim[, subs_pred[2]], na.rm = T),
                                length.out = lOut),
                         rep(mean(data_allEstim[, 'Pvalue'], na.rm = TRUE), lOut)), # rep(c(0, 1), each = lOut),
                         ncol = 2, byrow = FALSE)
  #Pred_data[, 4] <- Pred_data[,1]*Pred_data[,2]
  Pred_data <- as.data.frame(predict(mod_mv, newmods = Pred_data,  addx=TRUE))

  for(i in 1:length(names(Pred_data))){
    if(length(strsplit(names(Pred_data)[i], split = 'X.')[[1]]) == 1){
      names(Pred_data)[i] <- names(Pred_data)[i]
    } else {
      if(length(strsplit(names(Pred_data)[i], split = 'X.')[[1]]) == 2){
        tempName <- strsplit(names(Pred_data)[i], split = 'X.')[[1]][2]
        if(length(grep('Climate', tempName)) != 0){
          names(Pred_data)[i] <- strsplit(tempName, split = 'Temperature')[[1]]
        } else {
          names(Pred_data)[i] <- tempName
        }
      }
    }
  }

  if(byHemisphere){
    Pred_data %<>%
      dplyr::rename(Estimate = pred) %>%
      dplyr::mutate(., SD = (ci.ub - ci.lb) / 4,
                    Est_PlSD = Estimate + SD,
                    Est_MinSD = Estimate - SD,
                    Hemisphere = '')
  } else {
    Pred_data %<>%
      dplyr::rename(Estimate = pred) %>%
      dplyr::mutate(., SD = (ci.ub - ci.lb) / 4,
                    Est_PlSD = Estimate + SD,
                    Est_MinSD = Estimate - SD)
  }


  if(stats::coef(summary(mod_mv))$pval[2] < 0.05){  ## still hard-coded, decision on the line type
    lt = 1} else {lt = 2}
  if(byHemisphere){
    pl_CZ <- ggplot(data_allEstim, aes(x = data_allEstim[, names(Pred_data)[length(names(Pred_data))- 5]],
                                       y = Estimate, colour = Hemisphere)) +
      geom_point(alpha = 0.4) + theme_bw() +
      xlab(xLab) + ylab(yLab) +
      geom_line(data = Pred_data, aes(x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 5]],
                                      y = Estimate), linetype = lt, col = 'black') +
      geom_ribbon(data = Pred_data, aes(ymin = Est_MinSD, ymax = Est_PlSD,
                                        x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 5]]),
                  alpha=.2, col ='black') +
      scale_color_manual(values = c('deepskyblue2', 'goldenrod2')) +
      theme(legend.position = 'bottom') + ylim(miny, maxy)
  } else {
  pl_CZ <- ggplot(data_allEstim, aes(x = data_allEstim[, names(Pred_data)[length(names(Pred_data))- 4]],
                                     y = Estimate)) +
    geom_point(alpha = 0.4) + theme_bw() +
    xlab(xLab) + ylab(yLab) +
    geom_line(data = Pred_data, aes(x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 4]],
                                    y = Estimate), linetype = lt, col = 'black') +
    geom_ribbon(data = Pred_data, aes(ymin = Est_MinSD, ymax = Est_PlSD,
                                      x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 4]]),
                alpha=.2, col ='black') +
    theme(legend.position = 'bottom') + ylim(miny, maxy)
}
  if (!is.null(pdf_basename)) {
    grDevices::pdf(file = paste0(pdf_basename, '.pdf'))
  }
  print(pl_CZ)
  dev.off()

  return(list(Pred_data, pl_CZ))
}
