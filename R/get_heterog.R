#' Assess model heterogeneity
#'
#' \code{get_heterog} computes the heterogeneity metrics for meta-analysis.
#'
#' This function computes the heterogeneity metrics for the meta-analysis,
#' including I2, Q and p value associated with !
#'
#' @param mod A fitted mixed-effects meta-analytical model.
#' @param data Data frame used to fit the meta-analytical model,
#' which contains effect sizes and the respective sampling variances.
#'
#' @export
#' @return A tibble of the calculated heterogeneity metrics: Q reflects
#' the total amount of heterogeneity, Qp is the p value showing whether the observed
#' heterogeneity is larger than would be expected based on sampling variability,
#' I2 reflects the proportion of the total heterogeneity due to between-study
#' variance, ranges from 0 to 1, and I2_perRand is the proportion of the total
#' heterogeneity per random factor included in the model.

get_heterog <- function(mod, data){
  W <- diag(1/data$SError^2)
  X <- model.matrix(mod)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 <- sum(mod$sigma2) / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))
  I2_perRand <- mod$sigma2 / (sum(mod$sigma2) + (mod$k-mod$p)/sum(diag(P)))
  return(tibble::tibble(Q = mod$QE, Qp = mod$QEp,
                        I2 = I2, I2perRand = I2_perRand))
}
