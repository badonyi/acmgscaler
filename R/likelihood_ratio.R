# Main function and helpers for likelihood ratio calculation

density_ratio <- function(x) {
  # Computes positive likelihood ratios and their confidence interval

  # clamp to empirical range and rescale to [0,1]
  x <- clamp_rescale(x)

  # build a common density grid
  grid <- seq(0, 1, length.out = 1024)

  # map grid indices for each x
  idx <- findInterval(x, grid)

  # resample density grids for pathogenic and benign scores
  pd <- resample_density(x[names(x) == 'P'], grid)
  bd <- resample_density(x[names(x) == 'B'], grid)

  # compute variance-based penalty
  eps <- .Machine$double.eps
  mads <- apply(log(pd + eps) - log(bd + eps), 2, mad)
  lambda <- mads / max(mads) * sqrt(sum(diff(mads)^2)) / sum(mads)

  # regularise the log-LR ratio matrix
  llr_mat <- log(sweep(pd, 2, lambda, '+')) - log(sweep(bd, 2, lambda, '+'))

  # determine log-LR direction relative to x
  is_inv <- spearman(colMeans(llr_mat)[idx], x) < 0

  # monotonise the log-LR matrix in the principal direction
  for (i in seq_len(nrow(llr_mat))) {
    llr_mat[i, ] <- monotonise(llr_mat[i, ], is_inv)
  }

  # compute median and its confidence bounds
  llr_ci <- apply(llr_mat, 2, median_ci)

  # enforce shape constraint
  llr <- monotonise(llr_ci[2, ], is_inv)
  lo <- llr_ci[1, ] + (llr_ci[2, ] - llr)
  hi <- llr_ci[3, ] + (llr_ci[2, ] - llr)

  # exponentiate back and interpolate from the grid
  return(list(
    lower = exp(lo)[idx],
    lr = exp(llr)[idx],
    upper = exp(hi)[idx]
  ))
}


clamp_rescale <- function(x) {
  # Clamps to the range of 'P' and 'B' elements and rescales to [0,1]

  # limit to empirical range
  xpbr <- range(x[names(x) %in% c('P', 'B')], na.rm = TRUE, finite = TRUE)
  x[x < xpbr[1]] <- xpbr[1]
  x[x > xpbr[2]] <- xpbr[2]

  # min-max scale
  x <- (x - xpbr[1]) / (xpbr[2] - xpbr[1])

  return(x)
}


resample_density <- function(x, grid) {
  # Resamples densities by interpolating onto a predefined grid

  # determine bandwidth by biased cross-validation
  bw <- suppressWarnings(bw.bcv(x, nb = 1024))

  # preallocate a list for resamples
  boot_lst <- vector(mode = 'list', length = 1e3)
  n <- length(x)

  for (i in seq_len(1e3)) {
    tryCatch(
      expr = {
        d <- density(
          x = sample(x, n, replace = TRUE),
          n = 1024,
          bw = bw
        )

        # interpolate densities onto the common grid
        boot_lst[[i]] <- approx(d$x, d$y, xout = grid, rule = 2)$y
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
    expr = cor(x, y, method = 'spearman', use = 'pairwise.complete.obs'),
    warning = function(w) 0,
    error = function(e) 0
  )
}


monotonise <- function(x, flip) {
  # Bi-directional isotonic regression

  if (flip) {
    return(rev(isoreg(rev(x))$yf))
  } else {
    return(isoreg(x)$yf)
  }
}


median_ci <- function(x) {
  # Calculates the 95% confidence bounds for the population median
  # Conover, W. J. (John Wiley & Sons, 1999)

  n <- length(x)
  L <- qbinom(0.025, n, 0.5)
  U <- n - L + 1
  sx <- sort(x)
  c(lower = sx[L], median = median(x), upper = sx[n - L + 1])
}
