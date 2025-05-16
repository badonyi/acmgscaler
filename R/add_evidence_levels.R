# Classification helpers in add_likelihood_ratios()

add_evidence_levels <- function(df, pvst, lr_col) {
  # Helper for adding evidence classification used by all methods

  # derive evidence thresholds
  pp3_vst <- pvst
  pp3_st  <- pvst^0.5
  pp3_mo  <- pvst^0.25
  pp3_su  <- pvst^0.125
  bp4_su  <- 1 / pp3_su
  bp4_mo  <- 1 / pp3_mo
  bp4_st  <- 1 / pp3_st
  bp4_vst <- 1 / pvst

  # define evidence levels
  class_levels <- c(
    'b_very_strong',
    'b_strong',
    'b_moderate',
    'b_supporting',
    'indeterminate',
    'p_supporting',
    'p_moderate',
    'p_strong',
    'p_very_strong'
  )

  # preallocate a vector to store assigned classes
  out <- vector(mode = 'character', length = nrow(df))

  # assign classes based on LR confidence bounds
  lr_vec <- ifelse(
    test = df[[lr_col]] < 1,
    yes = df[[paste0(lr_col, '_upper')]],
    no = df[[paste0(lr_col, '_lower')]]
  )

  out[lr_vec <  bp4_vst] <- class_levels[1]
  out[lr_vec >= bp4_vst & lr_vec <  bp4_st]  <- class_levels[2]
  out[lr_vec >= bp4_st  & lr_vec <  bp4_mo]  <- class_levels[3]
  out[lr_vec >= bp4_mo  & lr_vec <  bp4_su]  <- class_levels[4]
  out[lr_vec >= bp4_su  & lr_vec <= pp3_su]  <- class_levels[5]
  out[lr_vec >  pp3_su  & lr_vec <= pp3_mo]  <- class_levels[6]
  out[lr_vec >  pp3_mo  & lr_vec <= pp3_st]  <- class_levels[7]
  out[lr_vec >  pp3_st  & lr_vec <= pp3_vst] <- class_levels[8]
  out[lr_vec >  pp3_vst] <- class_levels[9]
  out <- factor(out, levels = class_levels)

  # add to input data
  column <- paste0(gsub('_lr', '', lr_col), '_evidence')
  df[[column]] <- out

  return(df)
}


find_pvst <- function(prior) {
  # Helper to find the very strong threshold for pathogenicity given the prior

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
