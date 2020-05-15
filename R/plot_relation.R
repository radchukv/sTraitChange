#' Plot of the general relation (global effect size) and study-specific relations
#' (effect sizes)
#'
#' \code{plot_relation} Plots relations between independent and dependent variables
#'  displaying the effect sizes per each study and global effect across the studies
#'
#' @param fit_meta_data Tibble returned by the function \code{\link{fit_meta}}, that
#' contains in different data frames all the data necessary produce the plot.
#' @param rel Numeric specifying which relation (specified under the column 'names'
#' in the input tibble) to produce the plot for.
#' @param Covar Categorical specifying the name of the categorical variable that was included
#' as fixed-effect covariate in the meta-analysis. Defaults to NULL, in which case the
#' overall global effect size across all studies is plotted.
#' @param Demog_rate Character specifying the level of the demographic rate for which
#' analyses were conducted.
#' @param Trait_categ Character specifying the level of the trait, for which analyses were
#' conducted.
#' @param Clim Character specifying the level of the climatic variable, for which analyses
#' were conducted.
#' @param sel Character specifying the name to be included in the .pdf name, which indicates
#' the levels of demographic rate and trait for which the mixed-effect model was fitted.
#' @param colr Character vector specifying the colours to be used for the data points.
#' The length of the vector should correspond to the number of the levels in the
#' categorical explanatory variable included in the meta-analytical model, or should
#' be one (if the single global effect size across all studies is to be plotted).
#' @param folder_name Character specifyng the path to the directory in which
#' the results will be saved.
#' @param mar A vector specifying the plot margins, analogously to \code{\link[graphics]{par}}.
#'
#' @export
#'
#' @return Plots a relation between the dependent and independent variable
#' with a single line per effect size for each study, and global effect size(s)
#' together with their SE(s) shown with dashed line(s).
#'
#' @examples
#' Coefs_Aut <- readRDS(file = './output_forSEM_temp/PathCoefs_allMods_Temp_Weights_DD_Autocor.RDS')
#' meta_Phen_Repro <- fit_all_meta(data_MA = Coefs_Aut,
#'                                 Clim = 'Temperature',
#'                                 Demog_rate = 'Reproduction',
#'                                 Trait_categ = 'Phenological',
#'                                 tab = 'BirdType',
#'                                 Covar = NULL,
#'                                 sel = 'Phen_Repro',
#'                                 folder_name = './output_overall/',
#'                                 colr = c('black'))
#' plot_relation(fit_meta_data = meta_Phen_Repro,
#'                               rel = 1,
#'                               Covar = NULL,
#'                               Demog_rate = 'Reproduction',
#'                               Trait_categ = 'Phenological',
#'                               Clim = 'Temperature',
#'                               sel = 'Phen_Repro',
#'                               colr = c('black'),
#'                               folder_name  = './output_overall/',
#'                               mar = c(4, 4, 4, 2))
#'
#' meta_Phen_Repro_byTax <- fit_all_meta(data_MA = Coefs_Aut,
#'                                      Demog_rate = 'Reproduction',
#'                                      Trait_categ = 'Phenological',
#'                                      Clim = 'Temperature',
#'                                      tab = 'Taxon',
#'                                      Covar = 'Taxon',
#'                                      sel = 'Phen_Repro',
#'                                      folder_name = './output_overall/',
#'                                      colr = c('black', 'green', 'red'))
#' plot_relation(fit_meta_data = meta_Phen_Repro_byTax,
#'                               rel = x,
#'                               Covar = 'Taxon',
#'                               Demog_rate = 'Reproduction',
#'                               Trait_categ = 'Phenological',
#'                               Clim = 'Temperature',
#'                               sel = 'Phen_Repro',
#'                               colr = c('black', 'green', 'red'),
#'                               folder_name  = './output_overall/',
#'                               mar = c(4, 4, 7, 2))
#'
plot_relation <- function(fit_meta_data = meta_Phen_Repro_byTax,
                          rel = 1,
                          Covar = NULL,
                          Demog_rate = 'Survival',
                          Trait_categ = 'Phenological',
                          Clim = 'Temperature',
                          sel = 'Phen_surv',  ## selected combination of dem. rate and trait, for which the analyses are included in the data_globES
                          colr = c('black'),
                          folder_name = './output_overall/',
                          mar = c(3, 4, 7, 2)){

  ## extract datasets to work with
  data_ES <- fit_meta_data$data_meta[[1]]$data_EfS[[rel]]
  data_globES <- fit_meta_data$meta_res[[1]][rel, ]
  data_R2 <- fit_meta_data$data_meta[[1]]$data_R2[[rel]]
  if(! is.null(Covar)){
    data_CovarES <- fit_meta_data$data_meta[[1]]$data_Covar[[rel]]
  }
  ## start pdf if name if file defined
  coef <- plot_lab_name(Relation = data_globES[, 'Relation'],
                        Covar = Covar, Trait_categ = Trait_categ,
                        Clim = Clim, Demog_rate = Demog_rate)$pref_pdf
  if (!is.null(folder_name)) {
    grDevices::pdf(file = paste0(folder_name, sel, '_Slopes_', coef, '.pdf'))
  }

  graphics::par(mar = mar)

  plot(x = 0, y = 0, xlim = c(-1, 1), ylim = c(-1, 1), type = 'n',
       xlab = '', ylab = '')

  xax <- c(-1.1, 1.1) ## for the polygon

  for (i in 1:nrow(data_ES)){
    abline(a = 0, b = data_ES$Estimate[i], col = 'grey')
  }

  if(is.null(Covar)){
    if(as.numeric(data_globES[, 'pval_across']) < 0.05){
      LineT <- 1
    } else {
      LineT <- 2
    }
    yval <- 0 + as.numeric(data_globES[, 'Estimate'])*xax
    yvalU <- yval + as.numeric(data_globES[, 'SError'])
    yvalL <- yval - as.numeric(data_globES[, 'SError'])
    polygon(c(xax, rev(xax)), c(yvalL, rev(yvalU)), col = hsv(0.6,0.6,1, alpha = 0.4), border = FALSE)
    abline(a = 0, b = data_globES[, 'Estimate'], lwd = 6, col = colr, lty = LineT)

    graphics::mtext(paste0('Slope = ', round(data_globES[, 'Estimate'], 3)),
                    side = 3, line = 1, cex = 1.7, col = colr)
    if(as.numeric(data_globES[, 'pval_across']) < 0.0001){
      graphics::mtext(paste0('p < 0.0001'), side = 3, line = 0, cex = 1.7)
    }else{
      graphics::mtext(paste0('p = ', round(data_globES[, 'pval_across'], 4)), side = 3, line = 0, cex = 1.7)
    }
  } else {
    if(fit_meta_data$data_meta[[1]]$pval_Covar[rel] < 0.05){
      LineT <- 1
    } else {
      LineT <- 2
    }
    len_Covar <- length(unique(data_ES[, Covar]))
    for (i in 1:len_Covar){

      yval <- 0 + as.numeric(data_CovarES[i, 'Estimate'])*xax
      yvalU <- yval + as.numeric(data_CovarES[i, 'SError'])
      yvalL <- yval - as.numeric(data_CovarES[i, 'SError'])
      polygon(c(xax, rev(xax)), c(yvalL, rev(yvalU)),
              col = rgb(col2rgb(colr[i])[1], col2rgb(colr[i])[2], col2rgb(colr[i])[3], alpha = 0.6*255, maxColorValue = 255),
              border = FALSE)

      abline(a = 0, b = data_CovarES[i, 'Estimate'], lwd = 6, col = colr[i], lty = LineT)

      graphics::mtext(paste0(gsub(pattern = 'Continent', replacement = '', x= data_CovarES$Levels_Covar[i]),
                             ' slope = ', round(data_CovarES[i, 'Estimate'], 3)),
                      side = 3, line = 1 + i, cex = 1.7, col = colr[i])
    }
    legend('topright', legend = gsub(pattern = 'Continent', replacement = '', x= data_CovarES$Levels_Covar),
           col = colr, lwd= 3)
    if(fit_meta_data$data_meta[[1]]$pval_Covar[rel] < 0.0001){
      graphics::mtext(paste0('p < 0.0001'), side = 3, line = 0, cex = 1.7)
    }else{
      graphics::mtext(paste0('p = ', round(fit_meta_data$data_meta[[1]]$pval_Covar[rel], 4)), side = 3, line = 0, cex = 1.7)
    }
  }


  # Plot the value for Median R2
  if(regexpr(pattern = 'Tot', text = unique(data_R2$Relation)) < 0 &
     regexpr(pattern = 'Ind', text = unique(data_R2$Relation)) < 0){
    MedianR2 <- round(median(data_R2$R.squared, na.rm = TRUE), 2)
    graphics::text(x = -0.1, y = 0.7,
                   label = bquote(paste('Median R'^'2'*' = ', .(MedianR2))), cex = 1.7)
  }

  graphics::mtext(plot_lab_name(Relation = data_globES[, 'Relation'],
                                Covar = Covar, Trait_categ = Trait_categ,
                                Clim = Clim, Demog_rate = Demog_rate)$xlab_slope,
                  side = 1, line = 2, cex = 1.5)
  graphics::mtext(plot_lab_name(Relation = data_globES[, 'Relation'],
                                Covar = Covar, Trait_categ = Trait_categ,
                                Clim = Clim, Demog_rate = Demog_rate)$ylab_slope,
                  side = 2, line = 2, cex = 1.5)

  ## end pdf if name of the file specified
  if (!is.null(folder_name)) {
    grDevices::dev.off()
    message(paste0('a pdf named', paste0(folder_name, sel, '_Slopes_', coef, '.pdf'), ' has been created and saved!'))
  }
}
