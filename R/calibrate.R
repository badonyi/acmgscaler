#' Calibrate variant effect scores to ACMG/AMP evidence strength
#'
#'
#' @description
#' The function calculates the positive likelihood ratio (LR, equivalent to the
#' odds of pathogenicity) based on functional scores, e.g., from MAVEs or
#' computational predictors, and their truthset labels. Score intervals for
#' ACMG/AMP evidence levels are also computed. The input data requires at least
#' one numeric column with the score of interest and another column, named
#' \code{class}, with at least 10 pathogenic ('P') and 10 benign ('B') labels.
#' Different or missing labels are allowed, but will be renamed to 'U'.
#'
#'
#' @details
#' The function estimates the LR for each input score by resampling Gaussian
#' kernel density estimates of the pathogenic and benign score distributions.
#' Densities are mapped using linear interpolation and evaluated on a fixed-size
#' common grid. To stabilise the LRs in regions where densities approach zero, a
#' variance-based penalty is computed from log-LRs across 1,000 bootstrap
#' replicates. This penalty is used to regularise the log-LR matrix. The log-LRs
#' are monotonised in the principal direction of association with the input
#' scores. Final estimates for each score include the point estimate and its 95%
#' confidence interval. Score intervals for the different ACMG/AMP evidence
#' levels are interpolated from the grid based upon the confidence bounds.
#'
#'
#' @param df A dataframe. Must have a \code{class} column with values 'P'
#' (pathogenic) and 'B' (benign) labels, and a numeric column containing the
#' variant effect scores. At least 10 occurrences of each class are required.
#'
#' @param value (optional) A character string indicating the name of the numeric
#' column in \code{df} with the scores. If not provided, calibration will be
#' run on all numeric columns.
#'
#' @param prior A scalar in the range 0-1 representing the prior probability
#' of pathogenicity. Default 0.1.
#'
#' @param group (optional) A character string indicating the name of the column
#' with the grouping variable. Default NULL.
#'
#' @return A named list of dataframes. When grouping is not provided, the list
#' has a length of two where 'likelihood_ratios' is the input dataframe with
#' columns for LR and its confidence bounds (\code{column_name_lr},
#' \code{column_name_lr_lower} and \code{column_name_upper}). Assigned evidence
#' classifications can be found in the \code{evidence} column. The second
#' element in the list is named 'score_thresholds', which contain the lower and
#' upper bounds of the score interval for ACMG/AMP evidence levels. When
#' a grouping variable is provided, the returned object is a nested list with a
#' length equal to the unique group levels in the input data. Each of these
#' elements contain the 'likelihood_ratios' and 'score_thresholds' dataframes.
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
#' # load example data provided with the package
#' library(acmgscaler)
#' data(variant_data, package = 'acmgscaler')
#'
#' # compute likelihood ratios
#' calibrate(
#'   df = variant_data,
#'   # provide the column containing the variant scores (optional)
#'   value = 'score',
#'   # provide a grouping variable (optional)
#'   group = 'gene',
#'   # prior probability of pathogenicity (default 0.1)
#'   prior = 0.1
#' )
#'
#'
#' @export
calibrate <- function(df, value = NULL, prior = 0.1, group = NULL) {
  # S3 generic method for acmgscaler

  # set the seed locally
  user_seed <- get0('.Random.seed', envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (!is.null(user_seed)) {
      assign('.Random.seed', user_seed, envir = .GlobalEnv)
    } else {
      rm('.Random.seed', envir = .GlobalEnv)
    }
  }, add = TRUE)

  # set a new seed with R's current default
  set.seed(seed = 42, kind = 'Mersenne-Twister', sample.kind = 'Rejection')

  # check input suitability
  check_input(df = df, value = value, prior = prior, group = group)

  # method dispatch
  if (is.null(value)) {
    class(df) <- c('all_numeric', class(df))
  } else if (!is.null(group)) {
    class(df) <- c('grouped', class(df))
  } else {
    class(df) <- c('simple', class(df))
  }

  UseMethod('calibrate', df)
}


#' @inheritParams calibrate
#' @export
calibrate.simple <- function(df, value = NULL, prior = 0.1, group = NULL) {
  # Default method for caibrate()

  # find the very strong pathogenicity threshold for the prior
  pvst <- find_pvst(prior)

  # preserve input class
  class_in <- class(df)

  # label any NAs in class as unknown
  df$class <- toupper(as.character(df$class))
  df$class <- ifelse(!df$class %in% c('P', 'B'), 'U', df$class)

  # compute and add likelihood ratios
  lst_out <- build_grid(
    df = df,
    value = value,
    pvst = pvst
  )

  return(lst_out)
}


#' @inheritParams calibrate
#' @export
calibrate.all_numeric <- function(df, value = NULL, prior = 0.1, group = NULL) {
  # All-numeric-columns method for caibrate()

  # find the very strong pathogenicity threshold for the prior
  pvst <- find_pvst(prior)

  # label any NAs in class as unknown
  df$class <- toupper(as.character(df$class))
  df$class <- ifelse(!df$class %in% c('P', 'B'), 'U', df$class)

  # if grouped call grouped method
  if (!is.null(group)) {

    lst_out <- calibrate.grouped(
      df = df,
      value = value,
      prior = prior,
      group = group
    )

  } else {
    # otherwise run on all numeric columns
    lst_out <- build_grid_blockwise(
      df = df,
      value = value,
      pvst = pvst,
      group = group
    )
  }

  return(lst_out)
}


#' @inheritParams calibrate
#' @export
calibrate.grouped <- function(df, value = NULL, prior = 0.1, group = NULL) {
  # Grouped method for caibrate()

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
  lst_out <- lapply(seq_len(n), function(i) {

    cat('Group:', group_names[i], '|',  i, '/', n, '\n')

    build_grid_blockwise(
      df = split_dfs[[i]],
      value = value,
      pvst = pvst,
      group = group
    )
  })

  cat('\n')

  names(lst_out) <- group_names

  return(lst_out)
}


check_input <- function(df, value = NULL, prior = 0.1, group = NULL) {
  # Data check subroutine for calibrate()

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


build_grid_blockwise <- function(df, value, pvst, group) {
  # Wrapper for build_grid() for all_numeric and grouped methods

  # do not process if counts are insufficient within group
  class_n <- table(df$class)
  if (!all(c('P', 'B') %in% names(class_n)) || any(class_n[c('P', 'B')] < 10)) {

    message(paste(unique(df[[group]]), 'skipped due to low class counts.\n'))

    df[[paste0(value, '_lr_lower')]] <- NA_real_
    df[[paste0(value, '_lr')]] <- NA_real_
    df[[paste0(value, '_lr_upper')]] <- NA_real_
    df[[paste0(value, '_evidence')]] <- NA_real_

    out_lst <- list(
      likelihood_ratios = df,
      score_thresholds = build_threshold_df(value)
    )

    # match input class
    class(out_lst$score_thresholds) <- class(df)

    return(out_lst)
  }

  if (is.null(value)) {
    # process all numeric columns
    col_vec <- names(df)[sapply(df, is.numeric)]

    # preallocate a list for threshold tables
    threshold_lst <- vector(mode = 'list', length = length(col_vec))

    for (col in col_vec) {

      cat(sprintf('\r%-50s', paste(
        which(col == col_vec), '/', length(col_vec), '|', col
      )))

      # add likelihood ratios
      tmp_lst <- build_grid(
        df = df,
        value = col,
        pvst = pvst
      )

      # update df
      df <- tmp_lst$likelihood_ratios

      # add threshold table to the list
      threshold_lst[[which(col == col_vec)]] <- tmp_lst$score_thresholds
    }

    # join on evidence
    thresholds_merged <- Reduce(function(x, y) {
      merge(x, y, by = 'evidence', sort = FALSE)
    }, threshold_lst)

    # match input class
    class(thresholds_merged) <- class(df)

    # update output list
    out_lst <- list(
      likelihood_ratios = df,
      score_thresholds = thresholds_merged
    )

    cat('\n')

  } else {
    # process the user-defined column
    out_lst <- build_grid(
      df = df,
      value = value,
      pvst = pvst
    )
  }

  return(out_lst)
}


build_grid <- function(df, value, pvst) {
  # Wrapper for density_ratio() used by all methods

  # preserve input order
  df$original_order <- seq_len(nrow(df))

  # isolate rows with missing values
  is_na <- is.na(df[[value]])
  non_na_df <- df[!is_na, ]
  na_df <- df[is_na, ]

  # create output column names
  lr_lower_col <- paste0(value, '_lr_lower')
  lr_col <- paste0(value, '_lr')
  lr_upper_col <- paste0(value, '_lr_upper')
  evidence_col <- paste0(value, '_evidence')

  # do not process if not enough variants
  if (nrow(non_na_df) < 20) {
    df[[lr_lower_col]] <- NA_real_
    df[[lr_col]] <- NA_real_
    df[[lr_upper_col]] <- NA_real_
    df[[evidence_col]] <- NA_real_
    df$original_order <- NULL
    return(df)
  }

  # sort non-NA data by the value column
  sorted_df <- non_na_df[order(non_na_df[[value]], decreasing = TRUE), ]

  # derive evidence levels
  thresholds <- thresholds_from_pvst(pvst)

  # calculate likelihood ratios and score thresholds
  lr_lst <- density_ratio(
    x = stats::setNames(sorted_df[[value]], sorted_df$class),
    t = log(thresholds)
  )

  # add likelihood ratios to data
  sorted_df[[lr_lower_col]] <- lr_lst$lr_lower
  sorted_df[[lr_col]] <- lr_lst$lr
  sorted_df[[lr_upper_col]] <- lr_lst$lr_upper

  # classify based on likelihood ratio confidence interval
  sorted_df <- add_evidence_levels(
    df = sorted_df,
    pvst = unname(thresholds[8]),
    lr_col = lr_col
  )

  if (nrow(na_df) > 0) {
    na_df[[lr_lower_col]] <- NA_real_
    na_df[[lr_col]] <- NA_real_
    na_df[[lr_upper_col]] <- NA_real_
    na_df[[evidence_col]] <- NA_real_
  }

  # combine isolated dataframes
  df_out <- rbind(sorted_df, na_df)

  # arrange combined dataframe according to the original order
  df_out <- df_out[order(df_out$original_order), ]

  # remove the original order column
  df_out$original_order <- NULL

  # add score thresholds to evidence levels
  threshold_out <- build_threshold_df(value)
  threshold_out[[paste0(value, '_lower')]] <- lr_lst$t_lower
  threshold_out[[value]] <- lr_lst$t
  threshold_out[[paste0(value, '_upper')]] <- lr_lst$t_upper

  # match input class
  class(threshold_out) <- class(df)

  # return the annotated input data and the score intervals
  return(list(
    likelihood_ratios = df_out,
    score_thresholds = threshold_out)
  )
}


build_threshold_df <- function(value) {
  # Constructs an empty score_interval table

  threshold_df <- structure(
    list(
      evidence = factor(c(
        'Benign-VeryStrong|Benign-Strong',
        'Benign-Strong|Benign-Moderate',
        'Benign-Moderate|Benign-Supporting',
        'Benign-Supporting|Indeterminate',
        'Indeterminate|Pathogenic-Supporting',
        'Pathogenic-Supporting|Pathogenic-Moderate',
        'Pathogenic-Moderate|Pathogenic-Strong',
        'Pathogenic-Strong|Pathogenic-VeryStrong'
      ))
    ),
    class = 'data.frame',
    row.names = c(NA, -8L)
  )

  threshold_df[[paste0(value, '_lower')]] <- NA_real_
  threshold_df[[value]] <- NA_real_
  threshold_df[[paste0(value, '_upper')]] <- NA_real_

  return(threshold_df)
}
