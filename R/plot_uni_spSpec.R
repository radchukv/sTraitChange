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
#' @param yLab Character specifying the y axis label.
#' @param byHemisphere Boolean specifying whether raw data points should be coloured depending on
#' which Hemisphere they are coming from.
#' @param miny Numeric specifying the minimum limit for the y axis.
#' @param maxy Numeric specifying the maximum limit for the y axis.
#' @export
#' @importFrom magrittr "%>%"
#'
#' @return Plots a forest plot with effect sizes (and SEs) for each study,
#' and global effect size(s) on the bottom. The data used for the plot are
#' returned (invisibly).
#'
#' @examples
#' # prepare the dataset with phenological traits only
#' dataPaths_phen <- dataPaths %>%
#'     dplyr::filter(Trait_Categ == 'Phenological') %>%
#'     dplyr::mutate(absLat = abs(Latitude))
#' # fit the meta-analytical model
#' mod_CZ_PhenT_AbsLat <- metafor::rma.mv(Estimate ~ absLat + Pvalue,
#'                               V = Std.Error^2, random = list(~ 1|Species, ~1|ID, ~1|Location),
#'                               data = dataPaths_phen, method = 'ML')
#' plot_CZ_PhenT_AbsLat <- plot_uni_spSpec(data_allEstim = dataPaths_phen,
#'                                         mod_mv = mod_CZ_PhenT_AbsLat,
#'                                        lOut = 10, xLab = 'Absolute latitude',
#'                                        yLab = 'CZ estimate',
#'                                        pdf_basename = paste0(tempdir(), '/PlotCZ_PhenT_byAbsLat'),
#'                                        # attention: for this example we
#'                                        # write the data to a temporary directory,
#'                                        # to check its location type tempdir()
#'                                        byHemisphere = FALSE,
#'                                        miny = min(dataPaths_phen$Estimate) - 0.1,
#'                                        maxy = max(dataPaths_phen$Estimate) + 0.1)
#' message('Temporary directory is located at', tempdir())
#' message('Contents of the temporary directory after running fit_SEM()',
#' list.files(tempdir()))
#'
plot_uni_spSpec <- function(data_allEstim,
                            mod_mv,
                            lOut = 10,
                            pdf_basename = NULL,
                            xLab,
                            yLab = 'Estimate',
                            byHemisphere = FALSE,
                            miny = -1.2, maxy = 1.2){

  subs_pred <- names(as.data.frame(stats::predict(mod_mv,  addx=TRUE)))[seq(length(names(as.data.frame(stats::predict(mod_mv,  addx=TRUE)))) - nrow(stats::coef(summary(mod_mv))) + 1, length(names(as.data.frame(stats::predict(mod_mv,  addx=TRUE)))))]
  subs_pred <- unlist(lapply(strsplit(subs_pred, split = 'X.'), function(x){x[2]}))
  Pred_data <- matrix(c(seq(from = min(data_allEstim[, subs_pred[2]], na.rm = T),
                                to = max(data_allEstim[, subs_pred[2]], na.rm = T),
                                length.out = lOut),
                         rep(mean(data_allEstim[, 'Pvalue'], na.rm = TRUE), lOut)),
                         ncol = 2, byrow = FALSE)

  Pred_data <- as.data.frame(stats::predict(mod_mv, newmods = Pred_data,  addx=TRUE))

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
      dplyr::rename(Estimate = .data$pred) %>%
      dplyr::mutate(SD = (.data$ci.ub - .data$ci.lb) / 4,
                    Est_PlSD = .data$Estimate + .data$SD,
                    Est_MinSD = .data$Estimate - .data$SD,
                    Hemisphere = '')
  } else {
    Pred_data %<>%
      dplyr::rename(Estimate = .data$pred) %>%
      dplyr::mutate(SD = (.data$ci.ub - .data$ci.lb) / 4,
                    Est_PlSD = .data$Estimate + .data$SD,
                    Est_MinSD = .data$Estimate - .data$SD)
  }


  if(stats::coef(summary(mod_mv))$pval[2] < 0.01){  ## still hard-coded, decision on the line type
    ## the threshold is adjusted for multiple tests: 5 tests
    lt = 1 } else {lt = 2}
  Hemisphere <- Est_MinSD <- Est_PlSD <- NULL
  if(byHemisphere){
    Estimate <- NULL
    pl_CZ <- ggplot2::ggplot(data_allEstim, ggplot2::aes(x = data_allEstim[, names(Pred_data)[length(names(Pred_data))- 5]],
                                       y = Estimate, colour = Hemisphere)) +
      ggplot2::geom_point(alpha = 0.4) + ggplot2::theme_bw() +
      ggplot2::xlab(xLab) + ggplot2::ylab(yLab) +
      ggplot2::geom_line(data = Pred_data, ggplot2::aes(x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 5]],
                                      y = Estimate), linetype = lt, col = 'black') +
      ggplot2::geom_ribbon(data = Pred_data, ggplot2::aes(ymin = Est_MinSD, ymax = Est_PlSD,
                                        x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 5]]),
                  alpha=.2, col ='black') +
      ggplot2::scale_color_manual(values = c('deepskyblue2', 'goldenrod2')) +
      ggplot2::theme(legend.position = 'bottom',
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()) + ggplot2::ylim(miny, maxy)
  } else {
    Estimate <- NULL
  pl_CZ <- ggplot2::ggplot(data_allEstim, ggplot2::aes(x = data_allEstim[, names(Pred_data)[length(names(Pred_data))- 4]],
                                     y = Estimate)) +
    ggplot2::geom_point(alpha = 0.4) + ggplot2::theme_bw() +
    ggplot2::xlab(xLab) + ggplot2::ylab(yLab) +
    ggplot2::geom_line(data = Pred_data, ggplot2::aes(x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 4]],
                                    y = Estimate), linetype = lt, col = 'black') +
    ggplot2::geom_ribbon(data = Pred_data, ggplot2::aes(ymin = Est_MinSD, ymax = Est_PlSD,
                                      x = Pred_data[, names(Pred_data)[length(names(Pred_data))- 4]]),
                alpha=.2, col ='black') +
    ggplot2::theme(legend.position = 'bottom',
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::ylim(miny, maxy)
}
  if (!is.null(pdf_basename)) {
    path <- paste0(pdf_basename, '.pdf')
    grDevices::pdf(file = path)
  }
  print(pl_CZ)
  grDevices::dev.off()

  return(pl_CZ)
}
