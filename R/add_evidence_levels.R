# Classification helpers in calibrate()

add_evidence_levels <- function(df, pvst, lr_col) {
  # Helper for adding evidence levels in build_grid()

  # derive likelihood points for the given Pvst
  lp <- thresholds_from_pvst(pvst)

  # expand evidence with indeterminate
  evidence <- c(names(lp)[1:4], 'indeterminate', names(lp)[5:8])

  # assign evidence based on LR confidence bounds
  lr_vec <- ifelse(
    test = df[[lr_col]] < 1,
    yes = df[[paste0(lr_col, '_upper')]],
    no = df[[paste0(lr_col, '_lower')]]
  )

  # preallocate a vector to store assigned evidence
  out <- vector(mode = 'character', length = nrow(df))

  out[lr_vec <  lp[1]] <- evidence[1]
  out[lr_vec >= lp[1] & lr_vec <  lp[2]] <- evidence[2]
  out[lr_vec >= lp[2] & lr_vec <  lp[3]] <- evidence[3]
  out[lr_vec >= lp[3] & lr_vec <  lp[4]] <- evidence[4]
  out[lr_vec >= lp[4] & lr_vec <= lp[5]] <- evidence[5]
  out[lr_vec >  lp[5] & lr_vec <= lp[6]] <- evidence[6]
  out[lr_vec >  lp[6] & lr_vec <= lp[7]] <- evidence[7]
  out[lr_vec >  lp[7] & lr_vec <= lp[8]] <- evidence[8]
  out[lr_vec >  lp[8]] <- evidence[9]
  out <- factor(out, levels = evidence)

  # add to input data
  column <- paste0(gsub('_lr', '', lr_col), '_evidence')
  df[[column]] <- out

  return(df)
}


thresholds_from_pvst <- function(pvst) {
  # Helper to derive likelihood points for evidence levels given the Pvst

  c(
    'b_very_strong' = 1 / pvst,
    'b_strong' =      1 / pvst^0.5,
    'b_moderate' =    1 / pvst^0.25,
    'b_supporting' =  1 / pvst^0.125,
    'p_supporting' =  pvst^0.125,
    'p_moderate' =    pvst^0.25,
    'p_strong' =      pvst^0.5,
    'p_very_strong' = pvst
  )
}


thresholds_upon_prior <- function(prior) {
  # Helper to derive likelihood points for evidence levels given the prior

  thresholds_from_pvst(find_pvst(prior))
}


find_pvst <- function(prior) {
  # Helper to find the very strong threshold for pathogenicity for the prior

  # handle pvst bounds for prior extrema
  if (prior <= 0.01) return(8573)
  if (prior > 0.97) return(1)

  # preallocate vector of possible thresholds
  lr_vec <- vector(mode = 'numeric', length = 8573)
  for (pvst in 1:8573) {

    # derive supporting/moderate/strong
    su <- pvst^0.125
    mo <- pvst^0.25
    st <- pvst^0.5

    # likelihood ratios for clinical classes
    class_lr <- c(
      mo * pvst,             #1  LP_i
      mo * st,               #2  LP_ii
      su^2 * st,             #3  LP_iii
      mo^3,                  #4  LP_iv
      su^2 * mo^2,           #5  LP_v
      su^4 * mo,             #6  LP_vi
      st * pvst,             #7  P_ia
      mo^2 * pvst,           #8  P_ib
      su * mo * pvst,        #9  P_ic
      su^2 * pvst,           #10 P_id
      st^2,                  #11 P_ii
      mo^3 * st,             #12 P_iiia
      su^2 * mo^2 * st,      #13 P_iiib
      su^4 * mo * st         #14 P_iiic
    )

    # posterior probability of pathogenicity
    post_path <- (class_lr * prior) / ((class_lr - 1) * prior + 1)

    # criteria for LP/P
    lr_vec[pvst] <- sum(post_path[1:6] >= 0.9) + sum(post_path[7:14] >= 0.99)
  }

  # 13 out of 14 criteria met
  pvst <- which(lr_vec == 13)[1] - 1L

  return(pvst)
}


prettify_score_thresholds <- function(tbl) {
  # Colab-specific function for the interpretation of score intervals

  evidence_map <- c(
    'b_very_strong|b_strong'     = 'Benign-VeryStrong',
    'b_strong|b_moderate'        = 'Benign-Strong',
    'b_moderate|b_supporting'    = 'Benign-Moderate',
    'b_supporting|indeterminate' = 'Benign-Supporting',
    'indeterminate|p_supporting' = 'Pathogenic-Supporting',
    'p_supporting|p_moderate'    = 'Pathogenic-Moderate',
    'p_moderate|p_strong'        = 'Pathogenic-Strong',
    'p_strong|p_very_strong'     = 'Pathogenic-VeryStrong'
  )

  # rename evidence levels
  tbl$evidence <- evidence_map[match(tbl$evidence, names(evidence_map))]

  # find all valid columns
  col_check <- grep('_lower$', colnames(tbl), value = TRUE)
  col_ok <- col_check[sub('_lower$', '_upper', col_check) %in% names(tbl) &
                        sub('_lower$', '', col_check) %in% names(tbl)]

  pretty_threshold <- function(x, symbol) {
    # Helper to format thresholds

    is_valid <- !is.na(x)
    out <- rep(NA_character_, length(x))
    out[is_valid] <- sprintf('%#.5f', x[is_valid])
    symbol <- ifelse(is.na(symbol), '', symbol)
    out[is_valid] <- paste0(symbol[is_valid], out[is_valid])
    out
  }

  # iterate over all valid columns
  for (lower_col in col_ok) {
    upper_col <- sub('_lower$', '_upper', lower_col)
    score_col <- sub('_lower$', '', lower_col)

    lower_vals <- tbl[[lower_col]]
    upper_vals <- tbl[[upper_col]]

    decreasing <- is.unsorted(lower_vals, na.rm = TRUE)
    is_benign <- grepl('^Benign', tbl$evidence)

    threshold_vals <- rep(NA_real_, nrow(tbl))
    threshold_sym  <- rep(NA_character_, nrow(tbl))

    if (decreasing) {
      threshold_vals[!is_benign] <- lower_vals[!is_benign]
      threshold_sym[!is_benign]  <- '< '
      threshold_vals[is_benign]  <- upper_vals[is_benign]
      threshold_sym[is_benign]   <- '> '
    } else {
      threshold_vals[!is_benign] <- lower_vals[!is_benign]
      threshold_sym[!is_benign]  <- '> '
      threshold_vals[is_benign]  <- upper_vals[is_benign]
      threshold_sym[is_benign]   <- '< '
    }

    tbl[[score_col]] <- pretty_threshold(threshold_vals, threshold_sym)

    tbl[[lower_col]] <- NULL
    tbl[[upper_col]] <- NULL
  }

  colnames(tbl)[colnames(tbl) == 'evidence'] <- 'ACMG/AMP evidence strength'

  return(tbl)
}
