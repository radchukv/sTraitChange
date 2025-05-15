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
#' @param miny Numeric specifying the minimum limit for the y axis.
#' @param maxy Numeric specifying the maximum limit for the y axis.
#' @param col_var Character specifying what variable from the raw dataset (raw_dat)
#' is used to colour the lines.
#' @param lwd_leg Numeric specifying the line width in the legend.
#'
#' @inheritParams fit_all_meta
#' @inheritParams plot_forest
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#'
#' @return Plots a requested relation (CZ, ZG or CG) for each study in the dataset and
#' overlays the across-study global effect sizes estimated with the meta-analyses.
#' @examples
#' allYrs_T <- do.call('rbind', lapply(X = unique(dataSEM$ID), FUN = function(x){
#' subs <- droplevels(dataSEM[dataSEM$ID == x, ])
#' full_NA <- data.frame(Year = seq(min(subs$Year), max(subs$Year), by = 1), ID = x)}))
#'  consec_yrs_T <- merge(allYrs_T, dataSEM, by = c('ID','Year'), all= TRUE)
#'  # calculate GR
#'  temp_GR <- consec_yrs_T %>%
#'  dplyr::mutate(., Pop_mean_lag = c(Pop_mean[-1], NA)) %>%
#'  dplyr::mutate(., GR = log(Pop_mean_lag / Pop_mean)) %>%
#'  dplyr::filter(., !is.na(GR) & !is.na(Trait_mean) &
#'  !is.na(Demog_rate_mean) & !is.na(Pop_mean))
#' # get residuals for clim over time
#'  temp_GRRes <- split(temp_GR, temp_GR$ID) %>%
#'  purrr::map(., ~lm(Clim ~ Year, data = .)) %>%
#'  purrr::map2(.x = ., .y = split(temp_GR, f= temp_GR$ID),
#'  .f = ~broom::augment_columns(x = .x, data = .y)) %>%
#'  dplyr::bind_rows() %>%
#'  dplyr::select(-.rownames)
#'  # scale the variables as for the meta-analyses
#'  temp_std <- temp_GRRes %>%
#'  dplyr::group_by(ID) %>%
#'  dplyr::mutate(Trait_SE = Trait_SE / sd(Trait_mean, na.rm = TRUE),
#'                Demog_rate_SE = Demog_rate_SE /sd(Demog_rate_mean, na.rm = TRUE),
#'                det_Clim = as.numeric(scale(`.resid`)),
#'                Trait_mean = scale(Trait_mean),
#'                Demog_rate_mean = scale(Demog_rate_mean),
#'                Pop_mean = scale(Pop_mean),
#'                GR = scale(GR)) %>%
#'                dplyr::ungroup() %>%
#'                dplyr::mutate(Climatic_var = 'Temperature')
#'
#' temp_std$GR <- as.numeric(temp_std$GR[,1])
#'  # prepare the data with paths extracted
#'  allES_T <- dataPaths %>%
#'  dplyr::mutate(Climatic_var = 'Temperature',
#'  SError = Std.Error)
#'  wide_temp_all <- allES_T %>%
#'   dplyr::select(ID, Estimate, SError, P.Value, Relation) %>%
#'   tidyr::pivot_wider(id_cols = ID, names_from = Relation,
#'   values_from = c(Estimate, SError, P.Value), names_sep = '/')
#'
#'   metaD_temp_all <- allES_T %>%
#'   dplyr::distinct(ID, .keep_all = TRUE) %>%
#'   dplyr::select(-c(Estimate, SError, Relation,
#'                   P.Value, Pvalue, Count, Nyears,
#'                   WinDur, Ref.day, Ref.month,
#'                   WindowClose, deltaAIC, Trait_ageClass,
#'                   WeathQ, GenLength_y_IUCN))
#'
#' wide_tempES_all <- (merge(wide_temp_all, metaD_temp_all, by = 'ID'))
#'
#' # and estimating the global, across-study estimates
#' dataPaths_sp <- dataPaths %>%
#'                   dplyr::mutate(Species = dplyr::case_when(
#'                          Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
#'                          Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
#'                          Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
#'                          Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
#'                          Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
#'                          Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
#'                          TRUE ~ Species))
#'
#' dataPaths_sp$Species <- unlist(lapply(1:nrow(dataPaths_sp), FUN = function(x){
#'   binary <- strsplit(as.character(dataPaths_sp$Species[x]), " ")
#'   Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")}))
#' dataPaths_sp$Sp_phylo <- dataPaths_sp$Species
#'
#' meta_Phen_Cov <- fit_all_meta(data_MA = dataPaths_sp,
#'                               Demog_rate = NULL,
#'                               Trait_categ = 'Phenological',
#'                               Clim = 'Temperature',
#'                               Cov_fact = 'WeathQ',
#'                               COV = 'Pvalue',
#'                               sel = 'Temp_Phen_Cov',
#'                               folder_name = NULL,
#'                               colr = c('black', 'red'),
#'                               DD = 'n_effectGR',
#'                               simpleSEM = TRUE,
#'                               A = phyloMat,
#'                               all_Relations = c('Trait_mean<-det_Clim'))
#' globES_T <-  meta_Phen_Cov$meta_res[[1]] %>%
#' dplyr::filter(Levels_Covar == 'intrcpt') %>%
#' dplyr::mutate(Trait_Categ = 'Phenological',
#'             REL ='CZ')
#' # plot
#' PhenT_CZ <- plot_concept(Trait_categ = 'Phenological',
#'                         raw_dat = temp_std,
#'                         GlobES_dat = globES_T,
#'                         ES_dat = wide_tempES_all,
#'                         path = 'CZ',
#'                         xvar_raw = 'det_Clim',
#'                         yvar_raw = 'Trait_mean',
#'                         slope_ES = 'Estimate/Trait_mean<-det_Clim',
#'                         ylab = 'Phenology, Z',
#'                         xlab = 'Temperature, C',
#'                         miny = -4, maxy = 4)
plot_concept <- function(Trait_categ = 'Phenological',
                         raw_dat,
                         GlobES_dat,
                         ES_dat,
                         path = 'CZ',
                         xvar_raw = 'det_Clim',
                         yvar_raw = 'Trait_mean',
                         slope_ES = 'Estimate/Trait_mean<-det_Clim',
                         ylab = 'Trait', xlab = 'Climate',
                         miny = -6, maxy = 6,
                         col_var = NULL,
                         lwd_leg = 2){
  raw_dat <- raw_dat %>%
    dplyr::filter(.data$Trait_Categ == Trait_categ)
  GlobES_dat <- GlobES_dat[GlobES_dat$REL == path &
                             GlobES_dat$Trait_Categ == Trait_categ, ]
  GlobES_dat %<>%
    dplyr::mutate(ltype = dplyr::case_when(pval_Covar < 0.05 ~ '1',
                                           TRUE ~ '2'))
  ES_dat <- ES_dat %>%
    dplyr::filter(.data$Trait_Categ == Trait_categ)

    dat_rib <- data.frame(x = (seq(min(raw_dat$det_Clim),
                                        max(raw_dat$det_Clim),
                                        length.out = 10)))

    dat_rib %<>%
      dplyr::mutate(ymax = GlobES_dat$EfS_Upper * .data$x,
                    ymin = GlobES_dat$EfS_Low * .data$x,
                    Trait_mean = 0, GR = 0, Taxon = 'Bird')

  x <- ymin <- ymax <- Estimate <- ltype <- NULL
  if (requireNamespace("ggtext", quietly = TRUE)) {
    pl <- ggplot2::ggplot(raw_dat,
                          ggplot2::aes(x = .data[[xvar_raw]], y = .data[[yvar_raw]]),
                          alpha = 0.7) +
      ggplot2::lims(x = c(min(dat_rib$x), max(dat_rib$x)),
           y =  c(miny, maxy)) +
      ggplot2::geom_blank() +
      ggplot2::geom_abline(data = ES_dat,
                           ggplot2::aes(intercept = 0, slope = .data[[slope_ES]]),
                           alpha = 0.7, col = 'grey') +
      ggplot2::geom_ribbon(data = dat_rib,
                           ggplot2::aes(x = x, ymin= ymin, ymax = ymax),
                  fill = 'black', col = 'black',
                  alpha = 0.55) +
      ggplot2::geom_abline(data = GlobES_dat,
                           ggplot2::aes(intercept = 0, slope = Estimate,
                      lty = ltype), col = 'black', lwd = 1) +
      ggplot2::scale_colour_brewer(palette = 'Dark2') +
      ggplot2::scale_linetype_manual(values = c('1' = 1,
                                       '2' = 2),
                            labels = c('1' = 'p <= 0.1',
                                       '2' = 'p > 0.1')) +
      ggplot2::theme_bw() +
      ggplot2::ylab(ylab) + ggplot2::xlab(xlab) +
      ggplot2::theme(
            legend.position = 'bottom',  ## 22.06: to not interfere with other legends in the composite plot
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(size = 20),
            axis.text = ggplot2::element_text(size = 15),
            axis.title.x = ggtext::element_markdown(),
            axis.title.y = ggtext::element_markdown(),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 15),
            legend.key.width=unit(1,"cm"),
            plot.margin = margin(5.5,3,5.5,25, "pt")) +
      ggplot2::guides(lty = 'none')
  } else {
    message("to be able to produce this plot, you first must run install.packages('ggtext')!")
  }
  return(pl)
}
