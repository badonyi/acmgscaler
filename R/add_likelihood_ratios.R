#' Compute the likelihood ratio of pathogenicity
#'
#'
#' @description
#' This function calculates the likelihood ratio (equivalent to the odds of
#' pathogenicity) based on a numeric score, e.g., fitness scores from MAVE
#' experiments or computational predictors. The input data should have at least
#' one numeric column with the score of interest, and another column named
#' \code{class}, containing pathogenic ('P') and benign ('B') truthset labels
#' for the variants. Likelihood ratios will also be computed for unlabelled
#' variants or those different from 'P' or 'B' in the data.
#'
#'
#' @details
#' The function estimates the likelihood ratio of pathogenicity for each input
#' score by resampling Gaussian kernel density estimates of the 'P' (pathogenic)
#' and 'B' (benign) score distributions. Densities are evaluated on a common
#' grid and interpolated using linear interpolation. To stabilise the ratio in
#' regions where densities approach zero, a variance-based penalty is computed
#' from log-density ratios across 1,000 bootstrap replicates. This penalty is
#' used to regularise the log-likelihood ratio matrix. The log-likelihood ratios
#' are monotonised in the principal direction of association with the input
#' scores. Final estimates for each score include the point estimate (median)
#' and the 95% confidence interval for the location of the median, derived from
#' the resamples and mapped from the common grid to the input scores.
#'
#'
#' @param df A dataframe. Must contain a \code{class} column with values 'P'
#' (pathogenic) and 'B' (benign) variants, and a numeric column containing the
#' variant effect scores. At least 10 occurrences of each class is required.
#'
#' @param value (optional) A character string indicating the name of the numeric
#' column in \code{df} containing the scores. If not provided, likelihood ratios
#' will be computed for all numeric columns in the dataframe.
#'
#' @param prior A scalar in the range 0-1 representing the prior probability
#' of pathogenicity. Default 0.1.
#'
#' @param group (optional) A character string indicating the name of the column
#' with the grouping variable. Default NULL.
#'
#' @return The input dataframe with new columns for the likelihood ratio, in
#' the form \code{column_name_lr}, and \code{column_name_lr_lower} and
#' \code{column_name_upper} for the confidence bounds. The classifications
#' under the predefined prior are given in \code{column_name_evidence}.
#'
#'
#' @references Richards et al., 2015. Modeling the ACMG/AMP variant
#' classification guidelines as a Bayesian classification framework.
#' *Genetics in Medicine*.
#' [DOI: 10.1038/gim.2017.210](https://doi.org/10.1038/gim.2017.210)
#'
#' @references Tavtigian et al., 2018. Standards and guidelines for the
#' interpretation of sequence variants: a joint consensus recommendation of
#' the American College of Medical Genetics and Genomics and the Association
#' for Molecular Pathology.
#' *Genetics in Medicine*.
#' [DOI: 10.1038/gim.2015.30](https://doi.org/10.1038/gim.2015.30)
#'
#' @references Brnich et al., 2019. Recommendations for application of the
#' functional evidence PS3/BS3 criterion using the ACMG/AMP sequence variant
#' interpretation framework.
#' *Genome Medicine*.
#' [DOI: 10.1186/s13073-019-0690-2](https://doi.org/10.1186/s13073-019-0690-2)
#'
#' @references Pejaver et al., 2022. Calibration of computational tools for
#' missense variant pathogenicity classification and ClinGen recommendations
#' for PP3/BP4 criteria.
#' *The American Journal of Human Genetics*.
#' [DOI: 10.1016/j.ajhg.2022.10.013](https://doi.org/10.1016/j.ajhg.2022.10.013)
#'
#' @references van Loggerenberg et al., 2023. Systematically testing human HMBS
#' missense variants to reveal mechanism and pathogenic variation
#' *The American Journal of Human Genetics*.
#' [DOI: 10.1016/j.ajhg.2023.08.012](https://doi.org/10.1016/j.ajhg.2023.08.012)
#'
#'
#' @examples
#' # load the example data provided with the package
#' library(acmgscaler)
#' data(variant_data, package = 'acmgscaler')
#'
#' # compute likelihood ratios
#' add_likelihood_ratios(
#'   df = variant_data,
#'   # provide the column containing the variant scores (optional)
#'   value = 'score',
#'   # provide a grouping variable (optional)
#'   group = 'gene',
#'   # prior probability of pathogenicity
#'   prior = 0.1
#' )
#'
#'
#' @export
add_likelihood_ratios <- function(df,
                                  value = NULL,
                                  prior = 0.1,
                                  group = NULL) {
  # S3 generic method for likelihood ratio computation

  set.seed(seed = 42, kind = 'Mersenne-Twister', sample.kind = 'Rejection')

  check_input(
    df = df,
    value = value,
    prior = prior,
    group = group
  )

  # method dispatch
  if (is.null(value)) {
    class(df) <- c('all_numeric', class(df))
  } else if (!is.null(group)) {
    class(df) <- c('grouped', class(df))
  } else {
    class(df) <- c('simple', class(df))
  }

  UseMethod('add_likelihood_ratios', df)
}


#' @inheritParams add_likelihood_ratios
#' @export
add_likelihood_ratios.simple <- function(df,
                                         value = NULL,
                                         prior = 0.1,
                                         group = NULL) {
  # Default method for likelihood ratio computation

  # preserve input class
  class_in <- class(df)

  # label any NAs in class as unknown
  df$class <- toupper(as.character(df$class))
  df$class <- ifelse(!df$class %in% c('P', 'B'), 'U', df$class)

  # compute and add likelihood ratios
  df <- add_lr(
    df = df,
    value = value
  )

  # add resulting evidence levels
  df <- add_evidence_levels(
    df = df,
    pvst = find_pvst(prior),
    lr_col = paste0(value, '_lr')
  )

  rownames(df) <- NULL
  class(df) <- class_in
  return(df)
}


#' @inheritParams add_likelihood_ratios
#' @export
add_likelihood_ratios.all_numeric <- function(df,
                                              value = NULL,
                                              prior = 0.1,
                                              group = NULL) {
  # All-numeric-columns method for likelihood ratio computation

  # find very strong pathogenicity threshold for the prior
  pvst <- find_pvst(prior)

  # label any NAs in class as unknown
  df$class <- toupper(as.character(df$class))
  df$class <- ifelse(!df$class %in% c('P', 'B'), 'U', df$class)

  # if grouped, split and process by group
  if (!is.null(group)) {

    # split dataframe by group
    split_dfs <- split(df, df[[group]])
    n <- length(split_dfs)
    group_names <- names(split_dfs)

    # process each group individually
    classed_dfs <- lapply(seq_len(n), function(i) {

      cat('Group:', group_names[i], '|',  i, '/', n, '\n')

      add_lr_blockwise(
        df = split_dfs[[i]],
        value = value,
        pvst = pvst,
        group = group
      )
    })

    cat('\n')

    # bind lists into a dataframe
    out <- do.call(rbind, classed_dfs)

  } else {
    # if no group, run on all numeric columns
    out <- add_lr_blockwise(
      df = df,
      value = value,
      pvst = pvst,
      group = group
    )
  }

  if (is.null(out)) {
    return(df)
  }

  rownames(out) <- NULL
  class(out) <- class(df)
  return(out)
}


#' @inheritParams add_likelihood_ratios
#' @export
add_likelihood_ratios.grouped <- function(df,
                                          value = NULL,
                                          prior = 0.1,
                                          group = NULL) {
  # Grouped method for likelihood ratio computation

  # find the very strong pathogenicity threshold for the prior
  pvst <- find_pvst(prior)

  # label any NAs in class as unknown
  df$class <- toupper(as.character(df$class))
  df$class <- ifelse(!df$class %in% c('P', 'B'), 'U', df$class)

  # split dataframe by the grouping variable
  split_dfs <- split(df, df[[group]])
  n <- length(split_dfs)
  group_names <- names(split_dfs)

  # process each group individually
  classed_dfs <- lapply(seq_len(n), function(i) {

    cat('Group:', group_names[i], '|',  i, '/', n, '\n')

    add_lr_blockwise(
      df = split_dfs[[i]],
      value = value,
      pvst = pvst,
      group = group
    )
  })

  cat('\n')

  # bind lists into a dataframe
  out <- do.call(rbind, classed_dfs)

  rownames(out) <- NULL
  class(out) <- class(df)
  return(out)
}


check_input <- function(df,
                        value = NULL,
                        prior = 0.1,
                        group = NULL) {
  # Data check subroutine for add_likelihood_ratios()

  args <- mget(names(formals(sys.function())), envir = environment())
  if (!all(lengths(args[-1]) %in% c(0, 1))) {
    stop('All arguments except "df" must be scalar.')
  }

  if (!is.data.frame(df)) {
    stop('"df" must be a dataframe.')
  }

  if (!'class' %in% colnames(df)) {
    stop('"df" must contain a "class" column.')
  }

  # check if value exists
  if (!is.null(value) && !value %in% colnames(df)) {
    stop(paste0('"df" must contain a column named "', value, '".'))
  }

  if (!is.null(prior) && (!is.numeric(prior) || prior <= 0 || prior >= 1)) {
    stop('"prior" must be between 0 and 1.')
  }

  if (!is.null(group)) {
    if (!group %in% colnames(df)) {
      stop('Grouping variable ', group, ' not found in dataframe.')
    }

    if (anyNA(df[[group]])) {
      stop('The grouping variable "', group, '" contains missing values.')
    }
  }

  class_n <- table(df$class)
  if (!all(c('P', 'B') %in% names(class_n))) {
    stop('"class" column must contain both "P" and "B" levels.')
  } else if (any(class_n[c('P', 'B')] < 10)) {
    stop('Each class must contain at least 10 "P" and 10 "B" class instances.')
  }
}


add_lr_blockwise <- function(df, value, pvst, group) {
  # Wrapper for add_lr() for all.numeric and grouped methods

  # do not process if counts are insufficient
  class_n <- table(df$class)
  if (!all(c('P', 'B') %in% names(class_n)) || any(class_n[c('P', 'B')] < 10)) {

    message(paste0(unique(df[[group]]), '" skipped due to low class counts.\n'))

    df[[paste0(value, '_lr')]] <- NA_real_
    df[[paste0(value, '_lr_lower')]] <- NA_real_
    df[[paste0(value, '_lr_upper')]] <- NA_real_

    return(df)
  }

  if (is.null(value)) {
    # process all numeric columns
    col_vec <- names(df)[sapply(df, is.numeric)]

    for (col in col_vec) {

      cat(sprintf('\r%-50s', paste(
        which(col == col_vec), '/', length(col_vec), '|', col
      )))

      # add likelihood ratios
      df <- add_lr(
        df = df,
        value = col
      )

      # add resulting evidence levels
      df <- add_evidence_levels(
        df = df,
        pvst = pvst,
        lr_col = paste0(col, '_lr')
      )
    }

    cat('\n')

  } else {
    # process the user-defined column
    df <- add_lr(
      df = df,
      value = value
    )

    # add resulting evidence levels
    df <- add_evidence_levels(
      df = df,
      pvst = pvst,
      lr_col = paste0(value, '_lr')
    )
  }

  return(df)
}


add_lr <- function(df, value) {
  # Wrapper for density_ratio() used by all methods

  # preserve input order
  df$original_order <- seq_len(nrow(df))

  # isolate rows with missing values
  is_na <- is.na(df[[value]])
  non_na_df <- df[!is_na, ]
  na_df <- df[is_na, ]

  # create output column names
  lr_col <- paste0(value, '_lr')
  lr_lower_col <- paste0(value, '_lr_lower')
  lr_upper_col <- paste0(value, '_lr_upper')

  # do not process if not enough variants
  if (nrow(non_na_df) < 20) {
    df[[lr_col]] <- NA_real_
    df[[lr_lower_col]] <- NA_real_
    df[[lr_upper_col]] <- NA_real_
    df$original_order <- NULL
    return(df)
  }

  # sort non-NA data by the value column (faster in monotonise())
  sorted_df <- non_na_df[order(non_na_df[[value]], decreasing = TRUE), ]

  # calculate likelihood ratios
  lr_lst <- density_ratio(
    x = setNames(sorted_df[[value]], sorted_df$class)
  )

  # add likelihood ratios to data
  sorted_df[[lr_col]] <- lr_lst$lr
  sorted_df[[lr_lower_col]] <- lr_lst$lower
  sorted_df[[lr_upper_col]] <- lr_lst$upper

  if (nrow(na_df) > 0) {
    na_df[[lr_col]] <- NA_real_
    na_df[[lr_lower_col]] <- NA_real_
    na_df[[lr_upper_col]] <- NA_real_
  }

  # combine isolated dataframes
  df_out <- rbind(sorted_df, na_df)

  # arrange combined dataframe according to the original order
  df_out <- df_out[order(df_out$original_order), ]

  # remove the original order column
  df_out$original_order <- NULL

  return(df_out)
}
