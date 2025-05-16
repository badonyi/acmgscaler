#' Get score thresholds for ACMG/AMP evidence levels
#'
#'
#' @description
#' This function calculates classification thresholds for variant pathogenicity
#' based on their likelihood ratios. It returns the score values corresponding
#' to ACMG/AMP evidence strengths under the predefined prior.
#'
#'
#' @details
#' The function automatically detects columns ending in `_lr` and their
#' corresponding confidence bounds, and uses them to interpolate the score
#' thresholds. If a grouping variable is provided, the function computes
#' thresholds separately for each group and returns a list of dataframes.
#' The prior probability of pathogenicity can be specified, affecting the
#' positioning of evidence thresholds.
#'
#'
#' @param df A dataframe containing the variants to classify. Must contain a
#' score column (any name) and columns matching the `*_lr`, `*_lr_lower`,
#' and `*_lr_upper` patterns, added with \link{add_likelihood_ratios}.
#'
#' @param group (optional) A character string indicating the name of the column
#' with the grouping variable. Default NULL.
#'
#' @param prior (optional) A scalar in the range 0-1 representing the prior
#' probability of pathogenicity. Default 0.1.
#'
#' @return A dataframe (when group = NULL) or list of dataframes for groups
#' containing the thresholds for all score columns in the data. The first
#' column in the dataframe is `evidence`, which represents the threshold
#' between two adjacent evidence strengths, e.g., 'p_supporting|p_moderate'.
#'
#'
#' @examples
#' # load the example data provided with the package
#' library(acmgscaler)
#' data(variant_data, package = 'acmgscaler')
#'
#' # compute likelihood ratios
#' variant_data_lr <- add_likelihood_ratios(
#'   df = variant_data,
#'   # provide the column containing the variant scores (optional)
#'   value = 'score',
#'   # provide a grouping variable (optional)
#'   group = 'gene',
#'   # prior probability of pathogenicity for classification
#'   prior = 0.1
#' )
#'
#' # determine thresholds
#' get_thresholds(
#'   df = variant_data_lr,
#'   # provide the grouping variable
#'   group = 'gene',
#'   # prior probability of pathogenicity (should be the same as before)
#'   prior = 0.1
#' )
#'
#'
#' @export
get_thresholds <- function(df, group = NULL, prior = 0.1) {
  # Main function to derive score thresholds from data

  # check input
  check_get_thresholds(df, group = group, prior = prior)

  # detect likelihood ratio and score value columns
  lr_cols <- grep('_lr$', colnames(df), value = TRUE)
  value_cols <- gsub('_lr$', '', lr_cols)

  # ensure *_lr_lower and *_lr_upper columns exist for each score
  has_bounds <- function(x) {
    all(paste0(x, c('_lr_lower', '_lr_upper')) %in% colnames(df))
  }

  ok_value <- value_cols[vapply(value_cols, has_bounds, logical(1))]
  no_bounds <- setdiff(value_cols, ok_value)

  if (length(ok_value) < 1) {
    stop('No valid score variables with complete LR bounds found.')
  }

  if (length(no_bounds) > 0) {
    warning(
      'Missing LR bound columns for: ',
      paste(no_bounds, collapse = ', ')
    )
  }

  # keep only value columns with corresponding score values
  no_value <- ok_value[!ok_value %in% colnames(df)]
  value_cols <- setdiff(ok_value, no_value)

  if (length(value_cols) < 1) {
    stop(
      'Score variables missing: ',
      paste(ok_value, collapse = ', ')
    )
  }

  if (length(no_value) > 0) {
    warning(
      'Score variables missing: ',
      paste(no_value, collapse = ', ')
    )
  }

  # final LR column list based on validated value columns
  lr_cols <- paste0(value_cols, '_lr')

  if (!is.null(group)) {

    split_df <- split(df, df[[group]])

    return(lapply(split_df, function(df_group) {
      find_threshold(
        df = df_group,
        lr_cols = lr_cols,
        value_cols = value_cols,
        prior = prior
      )
    }))
  }

  return(find_threshold(
    df = df,
    lr_cols = lr_cols,
    value_cols = value_cols,
    prior = prior
  ))
}


check_get_thresholds <- function(df, group = NULL, prior = 0.1) {
  # Data check subroutine for get_thresholds()

  args <- mget(names(formals(sys.function())), envir = environment())
  if (!all(lengths(args[-1]) %in% c(0, 1))) {
    stop('All arguments except "df" must be scalar.')
  }

  if (!is.data.frame(df)) {
    stop('"df" must be a dataframe.')
  }

  if (!is.null(group)) {
    if (!group %in% colnames(df)) {
      stop('Grouping variable ', group, ' not found in dataframe.')
    }
  }

  triplet <- c('_lr$', '_lr_lower$', '_lr_upper$')
  if (!all(sapply(triplet, function(x) any(grepl(x, colnames(df)))))) {
    stop('Missing columns matching the pattern: *_lr, *_lr_lower, *_lr_upper.')
  }

  if (!is.null(prior) && (!is.numeric(prior) || prior <= 0 || prior >= 1)) {
    stop('"prior" must be between 0 and 1.')
  }
}


thresholds_upon_prior <- function(prior) {
  # Helper to derive evidence thresholds given a prior

  pvst <- find_pvst(prior)
  c(
    'b_very_strong' = 1 / pvst,
    'b_strong' = 1 / pvst^0.5,
    'b_moderate' = 1 / pvst^0.25,
    'b_supporting' = 1 / pvst^0.125,
    'p_supporting' = pvst^0.125,
    'p_moderate' = pvst^0.25,
    'p_strong' = pvst^0.5,
    'p_very_strong' = pvst
  )
}


safe_approx <- function(x, y, xout) {
  # Rule 1 approx() wrapper with error handling

  tryCatch(
    expr = suppressWarnings(approx(x, y, ties = mean, xout = xout, rule = 1)$y),
    error = function(e) rep(NA_real_, 8)
  )
}


find_threshold <- function(df, lr_cols, value_cols, prior) {
  # Wrapper for safe_approx() to construct the threshold table

  # allocate output dataframe
  threshold_out <- data.frame(
    evidence = c(
      'b_very_strong|b_strong',
      'b_strong|b_moderate',
      'b_moderate|b_supporting',
      'b_supporting|indeterminate',
      'indeterminate|p_supporting',
      'p_supporting|p_moderate',
      'p_moderate|p_strong',
      'p_strong|p_very_strong'
    )
  )

  # positive likelihood ratio thresholds for the given prior
  xout <- thresholds_upon_prior(prior)

  # iterate over all class columns
  for (i in seq_along(lr_cols)) {
    # interpolate thresholds
    score_col <- value_cols[i]
    score_lr <- ifelse(
      test = df[[lr_cols[i]]] < 1,
      yes = df[[paste0(lr_cols[i], '_upper')]],
      no = df[[paste0(lr_cols[i], '_lower')]]
    )

    threshold_out[[score_col]] <- safe_approx(
      x = score_lr,
      y = df[[score_col]],
      xout = xout
    )
  }

  # preserve input df class
  class(threshold_out) <- c(class(df), 'data.frame')

  return(threshold_out)
}
