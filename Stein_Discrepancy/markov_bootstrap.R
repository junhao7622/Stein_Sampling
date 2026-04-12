# markov_bootstrap.R

library(stats)  # Provides acf().

simulatepm <- function(N, p_change) {
  if (!is.numeric(N) || length(N) != 1 || N < 2) {
    stop("N must be a numeric scalar >= 2")
  }
  if (!is.numeric(p_change) || length(p_change) != 1 || p_change < 0 || p_change > 1) {
    stop("p_change must be a scalar in [0, 1]")
  }

  X          <- rep(-1, N)
  change_sign <- runif(N) < p_change
  for (i in seq_len(N)) {
    prev <- if (i == 1) X[N] else X[i - 1]
    if (change_sign[i]) {
      X[i] <- -prev
    } else {
      X[i] <-  prev
    }
  }
  return(X)
}

#' Estimate integrated autocorrelation time (IACT) from samples.
#'
#' Uses the initial-positive-sequence style truncation (sum positive ACF lags
#' until first non-positive value). This is used to set wild-bootstrap
#' Markov sign-flip probability for dependent chains.
#' @keywords internal
estimate_iact <- function(X, max_lag = NULL) {
  if (is.null(dim(X))) {
    X <- matrix(as.numeric(X), ncol = 1)
  } else {
    X <- as.matrix(X)
  }

  if (!is.numeric(X) || nrow(X) < 3) {
    return(list(iact = 1, iact_per_dimension = 1, max_lag = 0))
  }

  n <- nrow(X)
  d <- ncol(X)
  if (is.null(max_lag)) {
    max_lag <- min(n - 1, max(5, floor(sqrt(n))))
  }
  max_lag <- as.integer(max_lag)

  iact_dim <- rep(1, d)
  for (j in seq_len(d)) {
    xj <- X[, j]
    if (!all(is.finite(xj)) || stats::var(xj) <= 0) {
      iact_dim[j] <- 1
      next
    }

    acf_vals <- as.numeric(stats::acf(xj, lag.max = max_lag, plot = FALSE)$acf)[-1]
    if (length(acf_vals) == 0) {
      iact_dim[j] <- 1
      next
    }

    first_non_pos <- which(acf_vals <= 0)[1]
    if (is.na(first_non_pos)) {
      keep <- acf_vals
    } else if (first_non_pos == 1) {
      keep <- numeric(0)
    } else {
      keep <- acf_vals[seq_len(first_non_pos - 1)]
    }

    iact_dim[j] <- 1 + 2 * sum(keep)
    if (!is.finite(iact_dim[j]) || iact_dim[j] < 1) {
      iact_dim[j] <- 1
    }
  }

  list(
    iact = stats::median(iact_dim),
    iact_per_dimension = iact_dim,
    max_lag = max_lag
  )
}

#' Resolve Markov sign-flip probability.
#'
#' For dependent chains, Chwialkowski et al. (2016) recommend wild-bootstrap
#' block dependence aligned with sample dependence. We approximate this by
#' setting a_n ~= 1 / IACT when change_prob = "auto".
#' @keywords internal
resolve_markov_change_prob <- function(change_prob, X) {
  if (is.character(change_prob)) {
    if (length(change_prob) != 1 || !(change_prob %in% c("auto", "paper"))) {
      stop("change_prob must be numeric in [0, 1], 'auto', or 'paper'")
    }

    if (identical(change_prob, "paper")) {
      paper_info <- estimate_paper_change_prob(X)
      return(list(
        change_prob = paper_info$change_prob,
        auto = TRUE,
        mode = "paper",
        iact = paper_info$iact,
        iact_per_dimension = paper_info$iact_per_dimension,
        acf_lag_max = paper_info$acf_lag_max,
        thinning_suggested = paper_info$thinning_suggested,
        thinning_lag = paper_info$thinning_lag
      ))
    }

    iact_info <- estimate_iact(X)
    n <- if (is.null(dim(X))) length(X) else nrow(X)
    a_n <- max(1 / n, min(0.5, 1 / iact_info$iact))

    thinning_suggested <- is.finite(iact_info$iact) && iact_info$iact > 20
    if (thinning_suggested) {
      warning(
        paste0(
          "High chain autocorrelation detected (IACT ~= ",
          sprintf("%.2f", iact_info$iact),
          "). Consider thinning the chain before GOF testing."
        ),
        call. = FALSE
      )
    }

    return(list(
      change_prob = a_n,
      auto = TRUE,
        mode = "auto",
      iact = iact_info$iact,
      iact_per_dimension = iact_info$iact_per_dimension,
      acf_lag_max = iact_info$max_lag,
        thinning_suggested = thinning_suggested,
        thinning_lag = NA_integer_
    ))
  }

  if (!is.numeric(change_prob) || length(change_prob) != 1 || change_prob < 0 || change_prob > 1) {
    stop("change_prob must be numeric in [0, 1], 'auto', or 'paper'")
  }

  list(
    change_prob = as.numeric(change_prob),
    auto = FALSE,
    mode = "manual",
    iact = NA_real_,
    iact_per_dimension = NA_real_,
    acf_lag_max = NA_integer_,
    thinning_suggested = FALSE,
    thinning_lag = NA_integer_
  )
}

#' @keywords internal
estimate_paper_change_prob <- function(X) {
  if (is.null(dim(X))) {
    X <- matrix(as.numeric(X), ncol = 1)
  } else {
    X <- as.matrix(X)
  }

  n <- nrow(X)
  d <- ncol(X)
  lag_max <- min(n - 1, max(5, floor(sqrt(n))))

  thinning_lags <- rep(1L, d)
  for (j in seq_len(d)) {
    xj <- X[, j]
    if (!all(is.finite(xj)) || stats::var(xj) <= 0) {
      thinning_lags[j] <- 1L
      next
    }

    acf_vals <- as.numeric(stats::acf(xj, lag.max = lag_max, plot = FALSE)$acf)[-1]
    first_under_half <- which(abs(acf_vals) < 0.5)[1]
    thinning_lags[j] <- if (is.na(first_under_half)) lag_max else first_under_half
  }

  k <- max(1L, as.integer(stats::median(thinning_lags)))
  a_n <- max(1 / n, min(0.5, 0.1 / k))

  iact_info <- estimate_iact(X)
  thinning_suggested <- k > 1
  if (thinning_suggested) {
    warning(
      paste0(
        "Paper-style calibration selected: recommend thinning by k=",
        k,
        " to target lag-1 correlation below 0.5 before KSD testing."
      ),
      call. = FALSE
    )
  }

  list(
    change_prob = a_n,
    thinning_lag = k,
    acf_lag_max = lag_max,
    iact = iact_info$iact,
    iact_per_dimension = iact_info$iact_per_dimension,
    thinning_suggested = thinning_suggested
  )
}

# plot main
if (sys.nframe() == 0) {
  w <- simulatepm(1000, 0.02)
  print(acf(w, lag.max = 3, plot = FALSE)$acf)
  print(as.numeric(t(w) %*% w) / 10000)

  plot(w, type = "l", main = "simulatepm")
}
