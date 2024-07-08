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
                            yLab = 'Estimate'){

  subs_pred <- names(as.data.frame(predict(mod_mv,  addx=TRUE)))[seq(length(names(as.data.frame(predict(mod_mv,  addx=TRUE)))) -3, length(names(as.data.frame(predict(mod_mv,  addx=TRUE)))))]
  subs_pred <- unlist(lapply(strsplit(subs_pred, split = 'X.'), function(x){x[2]}))
  Pred_data <- matrix(c(rep(seq(from = min(data_allEstim[, subs_pred[1]], na.rm = T),
                                to = max(data_allEstim[, subs_pred[1]], na.rm = T),
                                length.out = lOut), 2),
                        rep(c(0, 1), each = lOut), rep(mean(data_allEstim[, 'Pvalue'], na.rm = TRUE), lOut*2),
                        rep(0, lOut*2)), ncol = 4, byrow = FALSE)
  Pred_data[, 4] <- Pred_data[,1]*Pred_data[,2]
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

  Pred_data %<>%
    dplyr::rename(Estimate = pred) %>%
    dplyr::mutate(., Climate = dplyr::recode(Climate,
                                             '0' = 'Precipitation', '1' = 'Temperature'),
                  SD = (ci.ub - ci.lb) / 4,
                  Est_PlSD = Estimate + SD,
                  Est_MinSD = Estimate - SD)


  pl_CZ <- ggplot(data_allEstim, aes(x = data_allEstim[, names(Pred_data)[length(names(Pred_data))- 6]],
                                     y = Estimate, group = Climate,
                                     col = Climate)) +
    geom_point(alpha = 0.4) + theme_bw() +
    xlab(xLab) + ylab(yLab) +
    geom_line(data = Pred_data, aes(x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 6]],
                                    y = Estimate, col = Climate)) +
    geom_ribbon(data = Pred_data, aes(ymin = Est_MinSD, ymax = Est_PlSD,
                                      x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 6]],
                                      fill = Climate), alpha=.2) +
    scale_color_manual(values = c('Temperature' ='red4', 'Precipitation' = 'royalblue4')) +
    scale_fill_manual(values = c('Temperature' = 'tomato1', 'Precipitation' = 'lightslateblue')) +
    theme(legend.position = 'bottom')

  if (!is.null(pdf_basename)) {
    grDevices::pdf(file = paste0(pdf_basename, '.pdf'))
  }
  print(pl_CZ)
  dev.off()

  return(list(Pred_data, pl_CZ))
}
