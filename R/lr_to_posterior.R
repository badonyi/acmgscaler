#' Conversion between likelihood ratio and posterior probability
#'
#'
#'#' @description
#' This function converts positive likelihood ratios into posterior probability
#' of pathogenicity.
#'
#'
#' @param lr A numeric vector of likelihood ratios.
#'
#' @param prior A scalar in the range 0-1 representing the prior probability
#' of pathogenicity. Default 0.1.
#'
#' @return A numeric vector of the same length as `lr` containing the posterior
#' probability of pathogenicity values.
#'
#'
#' @examples
#' library(acmgscaler)
#' lr_to_posterior(lr = c(3.5, 35, 350), prior = 0.1)
#'
#'
#' @export
lr_to_posterior <- function(lr, prior = 0.1) {
  posterior_odds <- (prior / (1 - prior)) * lr
  posterior_odds / (1 + posterior_odds)
}
