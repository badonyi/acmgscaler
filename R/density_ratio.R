# Main function and helpers for density ratio calculation

density_ratio <- function(x, t) {
  # Computes positive likelihood ratios and their confidence interval

  # empirical range of x
  xpbr <- range(x[names(x) %in% c('P', 'B')], na.rm = TRUE, finite = TRUE)

  # clamp to empirical range and rescale to [0,1]
  x <- clamp_rescale(x, xpbr)

  # build a common density grid
  grid <- seq(0, 1, length.out = 1024)

  # inverse map to original scale
  grid_x <- grid * (xpbr[2] - xpbr[1]) + xpbr[1]

  # map grid indices for each x
  idx <- findInterval(x, grid)

  # resample density grids for pathogenic and benign scores
  pd <- resample_density(x[names(x) == 'P'], grid)
  bd <- resample_density(x[names(x) == 'B'], grid)

  # compute variance-based penalty
  eps <- .Machine$double.eps
  mads <- apply(log(pd + eps) - log(bd + eps), 2, mad)
  lambda <- mads / max(mads) * sqrt(sum(diff(mads)^2)) / sum(mads)

  # regularise the log-LR matrix
  llr_mat <- log(sweep(pd, 2, lambda, '+')) - log(sweep(bd, 2, lambda, '+'))

  # determine log-LR direction relative to x
  is_inv <- spearman(colMeans(llr_mat)[idx], x) < 0

  # monotonise and smooth the log-LR matrix in the principal direction
  for (i in seq_len(nrow(llr_mat))) {
    llr_mat[i, ] <- monotonise(llr_mat[i, ], is_inv)
  }

  # compute 90% percentile-based confidence intervals
  llr_ci <- apply(llr_mat, 2, median_ci)

  # strictly enforce shape constraint on estimates
  llr <- monotonise(llr_ci[2, ], is_inv)
  lo <- llr_ci[1, ] + llr_ci[2, ] - llr
  hi <- llr_ci[3, ] + llr_ci[2, ] - llr

  # map score intervals from the grid back onto x
  x_t <- score_interval(lo, llr, hi, grid_x, is_inv, t)

  # exponentiate log-LRs and map x indices on the grid
  return(list(
    lr_lower = exp(lo)[idx],
    lr = exp(llr)[idx],
    lr_upper = exp(hi)[idx],
    t_lower = x_t$lower,
    t = x_t$x,
    t_upper = x_t$upper
  ))
}


clamp_rescale <- function(x, xpbr) {
  # Clamps to the range of 'P' and 'B' elements and rescales to [0,1]

  # limit to empirical range
  x[x < xpbr[1]] <- xpbr[1]
  x[x > xpbr[2]] <- xpbr[2]

  # min-max scale
  x <- (x - xpbr[1]) / (xpbr[2] - xpbr[1])

  return(x)
}


resample_density <- function(x, grid) {
  # Resamples densities by interpolating onto a predefined grid

  # determine bandwidth by biased cross-validation
  bw <- suppressWarnings(stats::bw.bcv(x, nb = 1024))

  # preallocate a list for resamples
  boot_lst <- vector(mode = 'list', length = 1e3)
  n <- length(x)

  for (i in seq_len(1e3)) {
    tryCatch(
      expr = {
        d <- stats::density(
          x = sample(x, n, replace = TRUE),
          n = 1024,
          bw = bw
        )

        # interpolate densities onto the common grid
        boot_lst[[i]] <- stats::approx(d$x, d$y, xout = grid, rule = 2)$y
      },

      error = function(e) {
        boot_lst[[i]] <- NULL
      })
  }

  # return a matrix of density resamples
  return(do.call(rbind, boot_lst))
}


spearman <- function(x, y) {
  # Spearman correlation wrapper

  tryCatch(
    expr = stats::cor(x, y, method = 'spearman', use = 'pairwise.complete.obs'),
    warning = function(w) 0,
    error = function(e) 0
  )
}


monotonise <- function(x, is_inv) {
  # Bi-directional isotonic regression

  if (is_inv) {
    return(rev(stats::isoreg(rev(x))$yf))
  } else {
    return(stats::isoreg(x)$yf)
  }
}


median_ci <- function(x) {
  # Calculates the 95% confidence bounds for the population median
  # Conover, W. J. (John Wiley & Sons, 1999)

  n <- length(x)
  L <- stats::qbinom(0.025, n, 0.5)
  U <- n - L + 1
  sx <- sort(x)
  c(lower = sx[L], median = stats::median(x), upper = sx[n - L + 1])
}


lr_to_x <- function(x, y, xout) {
  # Rule 1 stats::approx() wrapper to determine score intervals

  tryCatch(
    expr = suppressWarnings(
      stats::approx(x, y, ties = mean, xout = xout, rule = 1)$y
    ),
    error = function(e) rep(NA_real_, 8)
  )
}


score_interval <- function(lo, llr, hi, grid_x, is_inv, t) {
  # Map score intervals from the LR gird back onto x

  # score thresholds on original scale
  x_lo <- lr_to_x(lo, grid_x, t)
  x_th <- lr_to_x(llr, grid_x, t)
  x_hi <- lr_to_x(hi, grid_x, t)
  names(x_lo) <- names(x_th) <- names(x_hi) <- names(t)

  return(list(
    lower = if (is_inv) x_lo else x_hi,
    x = x_th,
    upper = if (is_inv) x_hi else x_lo
  ))
}
