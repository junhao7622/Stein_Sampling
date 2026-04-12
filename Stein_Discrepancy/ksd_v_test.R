# ksd_v_test.R

if (file.exists("stein_helpers.R")) {
  source("stein_helpers.R")
}
if (file.exists("markov_bootstrap.R")) {
  source("markov_bootstrap.R")
}

#' Kernel Stein Discrepancy GOF Test (V-statistics)
#'
#' Performs a goodness-of-fit test using a Gaussian-kernel Stein discrepancy
#' with V-statistics. Supports an aggregated global test or per-dimension tests
#' with Holm correction.
#'
#' @param X Numeric vector or matrix of samples. A vector is treated as `n x 1`.
#' @param score_function Function of `X` returning score values with shape `n x d`
#'   (or a length-`n` vector when `d = 1`).
#' @param boot_method One of `"rademacher"` or `"markov"`.
#' @param scaling Positive kernel scaling parameter. If `NULL`, uses the median
#'   squared-distance heuristic.
#' @param nboot Number of bootstrap samples.
#' @param median_max_samples Maximum sample size used by median heuristic.
#' @param median_use_sampling If `TRUE`, median heuristic uses subsampling above
#'   `median_max_samples`.
#' @param median_seed Random seed used by heuristic sampling.
#' @param block_size Optional block size for matrix-free large-N evaluation.
#' @param block_threshold Auto-enables block mode for `n > block_threshold`.
#' @param test_type One of `"aggregated"` or `"per_dimension"`.
#' @param change_prob Markov sign-flip probability when `boot_method = "markov"`.
#'   Use `"auto"` to set `a_n ~= 1 / IACT` from the sample chain.
#'   For `test_type = "per_dimension"`, scaling is re-estimated for each
#'   1D slice and scalar `scaling` input is ignored.
#' @param return_matrices Logical flag. If `TRUE`, stores per-dimension
#'   `n x n` matrices in output info (high memory cost).
#'
#' @return An object of class `htest` with test statistic(s), p-value information,
#' bootstrap samples, and diagnostic metadata in `info`.
#'
#' @export
ksd_v_test <- function(X, score_function,
                       boot_method = c("rademacher", "markov"),
                       scaling = NULL,
                       nboot = 1000,
                       median_max_samples = 2000,
                       median_use_sampling = TRUE,
                       median_seed = 123,
                       block_size = NULL,
                       block_threshold = 5000,
                       test_type = c("aggregated", "per_dimension"),
                       change_prob = "auto",
                       return_matrices = FALSE) {
  data_name <- deparse(substitute(X))
  test_type <- match.arg(test_type)
  boot_method <- match.arg(boot_method)
  scaling_input <- scaling

  x_mat <- as_samples_matrix(X)
  n <- nrow(x_mat)
  d <- ncol(x_mat)

  if (!is.numeric(nboot) || length(nboot) != 1 || nboot <= 0) {
    stop("nboot must be a positive scalar")
  }
  nboot <- as.integer(nboot)

  if (!is.null(block_size)) {
    if (!is.numeric(block_size) || length(block_size) != 1 || block_size < 2) {
      stop("block_size must be NULL or a numeric scalar >= 2")
    }
    block_size <- as.integer(block_size)
  }

  if (!is.numeric(block_threshold) || length(block_threshold) != 1 || block_threshold < 2) {
    stop("block_threshold must be a numeric scalar >= 2")
  }
  block_threshold <- as.integer(block_threshold)

  grads <- as_gradient_matrix(score_function, x_mat)

  if (is.null(scaling)) {
    scaling <- find_median_distance(
      x_mat,
      max_samples = median_max_samples,
      use_sampling = median_use_sampling,
      seed = median_seed
    )
  }

  if (!is.numeric(scaling) || length(scaling) != 1 || !is.finite(scaling) || scaling <= 0) {
    stop("scaling must be a positive scalar")
  }

  cp_info <- resolve_change_probability_for_bootstrap(boot_method, change_prob, x_mat)
  block_info <- resolve_block_settings(n, block_size, block_threshold)

  if (test_type == "aggregated") {
    if (block_info$use_block_mode) {
      blocked <- compute_v_stat_bootstrap_blocked(
        X = x_mat,
        grads = grads,
        scaling = scaling,
        boot_method = boot_method,
        nboot = nboot,
        change_prob = cp_info$change_prob,
        block_size = block_info$block_size,
        dim_index = NULL
      )

      result <- list(
        statistic = c(ksd_v = blocked$stat_value),
        p.value = blocked$p_value,
        method = paste("Kernel Stein Discrepancy (V-statistics) -", test_type),
        data.name = data_name,
        bootstrap_samples = blocked$boot_stats,
        info = list(
          scaling = scaling,
          scaling_per_dimension = scaling,
          U_matrix = NULL,
          nboot = nboot,
          boot_method = boot_method,
          change_prob = cp_info$change_prob,
          auto_change_prob = cp_info$auto,
          change_prob_mode = if (!is.null(cp_info$mode)) cp_info$mode else "manual",
          iact = cp_info$iact,
          iact_per_dimension = cp_info$iact_per_dimension,
          thinning_suggested = cp_info$thinning_suggested,
          thinning_lag = if (!is.null(cp_info$thinning_lag)) cp_info$thinning_lag else NA_integer_,
          block_mode = TRUE,
          block_size = block_info$block_size,
          max_offdiag_kernel = blocked$max_offdiag_kernel
        )
      )
    } else {
      u_mat <- compute_v_matrix(x_mat, grads, scaling)
      boot <- perform_bootstrap(u_mat, boot_method, nboot, cp_info$change_prob)

      result <- list(
        statistic = c(ksd_v = n * mean(u_mat)),
        p.value = boot$p_value,
        method = paste("Kernel Stein Discrepancy (V-statistics) -", test_type),
        data.name = data_name,
        bootstrap_samples = boot$boot_stats,
        info = list(
          scaling = scaling,
          scaling_per_dimension = scaling,
          U_matrix = u_mat,
          nboot = nboot,
          boot_method = boot_method,
          change_prob = cp_info$change_prob,
          auto_change_prob = cp_info$auto,
          change_prob_mode = if (!is.null(cp_info$mode)) cp_info$mode else "manual",
          iact = cp_info$iact,
          iact_per_dimension = cp_info$iact_per_dimension,
          thinning_suggested = cp_info$thinning_suggested,
          thinning_lag = if (!is.null(cp_info$thinning_lag)) cp_info$thinning_lag else NA_integer_,
          block_mode = FALSE,
          block_size = NA_integer_,
          max_offdiag_kernel = NA_real_
        )
      )
    }
  } else {
    if (!is.null(scaling_input)) {
      warning(
        "Input 'scaling' is ignored for test_type = 'per_dimension'; using per-dimension median heuristic.",
        call. = FALSE
      )
    }

    stats <- numeric(d)
    p_vals <- numeric(d)
    scaling_dim <- numeric(d)
    boot_mat <- matrix(NA_real_, nrow = nboot, ncol = d)
    if (!is.logical(return_matrices) || length(return_matrices) != 1 || is.na(return_matrices)) {
      stop("return_matrices must be TRUE or FALSE")
    }
    if (return_matrices && block_info$use_block_mode) {
      warning(
        "return_matrices=TRUE is not available in block mode; returning NULL matrices to avoid memory explosion.",
        call. = FALSE
      )
    }
    u_mats <- if (return_matrices && !block_info$use_block_mode) vector("list", d) else NULL
    max_offdiag_kernel <- numeric(d)

    for (dim_index in seq_len(d)) {
      x_dim_mat <- matrix(x_mat[, dim_index], ncol = 1)
      scaling_dim[dim_index] <- find_median_distance(
        x_dim_mat,
        max_samples = median_max_samples,
        use_sampling = median_use_sampling,
        seed = if (is.null(median_seed)) NULL else as.integer(median_seed) + dim_index - 1
      )
      if (!is.numeric(scaling_dim[dim_index]) || !is.finite(scaling_dim[dim_index]) || scaling_dim[dim_index] <= 0) {
        stop("Per-dimension scaling must be positive")
      }

      if (block_info$use_block_mode) {
        blocked <- compute_v_stat_bootstrap_blocked(
          X = x_mat,
          grads = grads,
          scaling = scaling_dim[dim_index],
          boot_method = boot_method,
          nboot = nboot,
          change_prob = cp_info$change_prob,
          block_size = block_info$block_size,
          dim_index = dim_index
        )
        stats[dim_index] <- blocked$stat_value
        p_vals[dim_index] <- blocked$p_value
        boot_mat[, dim_index] <- blocked$boot_stats
        if (!is.null(u_mats)) {
          u_mats[[dim_index]] <- NULL
        }
        max_offdiag_kernel[dim_index] <- blocked$max_offdiag_kernel
      } else {
        u_dim <- compute_v_matrix_single_dimension(
          x_mat,
          grads,
          scaling_dim[dim_index],
          dim_index,
          k_mat = NULL
        )
        stats[dim_index] <- n * mean(u_dim)
        boot <- perform_bootstrap(u_dim, boot_method, nboot, cp_info$change_prob)
        p_vals[dim_index] <- boot$p_value
        boot_mat[, dim_index] <- boot$boot_stats
        if (!is.null(u_mats)) {
          u_mats[[dim_index]] <- u_dim
        }
        max_offdiag_kernel[dim_index] <- max_offdiag_from_matrix(rbf_kernel_matrix(x_dim_mat, scaling_dim[dim_index]))
        if (is.null(u_mats)) {
          rm(u_dim)
        }
      }
    }

    names(stats) <- paste0("dim_", seq_len(d))
    names(p_vals) <- paste0("dim_", seq_len(d))
    p_vals_corrected <- p.adjust(p_vals, method = "holm")

    result <- list(
      statistic = stats,
      p.values_corrected = p_vals_corrected,
      method = paste("Kernel Stein Discrepancy (V-statistics) -", test_type),
      data.name = data_name,
      bootstrap_samples = boot_mat,
      info = list(
        scaling = scaling_dim,
        scaling_per_dimension = scaling_dim,
        U_matrices = u_mats,
        return_matrices = !is.null(u_mats),
        raw_p_values = p_vals,
        nboot = nboot,
        boot_method = boot_method,
        change_prob = cp_info$change_prob,
        auto_change_prob = cp_info$auto,
        change_prob_mode = if (!is.null(cp_info$mode)) cp_info$mode else "manual",
        iact = cp_info$iact,
        iact_per_dimension = cp_info$iact_per_dimension,
        thinning_suggested = cp_info$thinning_suggested,
        thinning_lag = if (!is.null(cp_info$thinning_lag)) cp_info$thinning_lag else NA_integer_,
        block_mode = block_info$use_block_mode,
        block_size = if (block_info$use_block_mode) block_info$block_size else NA_integer_,
        max_offdiag_kernel = max_offdiag_kernel
      )
    )
  }

  class(result) <- "htest"
  result
}

#' @keywords internal
resolve_change_probability_for_bootstrap <- function(boot_method, change_prob, x_mat) {
  if (boot_method != "markov") {
    return(list(
      change_prob = NA_real_,
      auto = FALSE,
      iact = NA_real_,
      iact_per_dimension = NA_real_,
      thinning_suggested = FALSE
    ))
  }

  resolve_markov_change_prob(change_prob, x_mat)
}

#' @keywords internal
resolve_block_settings <- function(n, block_size, block_threshold) {
  auto_block <- n > block_threshold
  if (is.null(block_size)) {
    block_size <- if (auto_block) min(1024L, n) else n
  }

  list(
    use_block_mode = auto_block || block_size < n,
    block_size = as.integer(block_size)
  )
}

#' @keywords internal
as_samples_matrix <- function(x) {
  if (is.null(dim(x))) {
    x <- matrix(as.numeric(x), ncol = 1)
  } else {
    x <- as.matrix(x)
  }

  if (!is.numeric(x)) {
    stop("X must be numeric")
  }
  if (nrow(x) < 2) {
    stop("X must contain at least two samples")
  }

  x
}

#' @keywords internal
as_gradient_matrix <- function(score_function, x_mat) {
  n <- nrow(x_mat)
  d <- ncol(x_mat)

  grads <- score_function(x_mat)

  if (is.null(dim(grads))) {
    if (d != 1 || length(grads) != n) {
      stop("score_function output must be n x d (or length n when d = 1)")
    }
    grads <- matrix(as.numeric(grads), ncol = 1)
  } else {
    grads <- as.matrix(grads)
  }

  if (nrow(grads) != n || ncol(grads) != d) {
    stop("score_function output shape must match X (n x d)")
  }

  grads
}

#' @keywords internal
squared_distance_matrix <- function(x_mat) {
  x_norm <- rowSums(x_mat * x_mat)
  x2 <- matrix(x_norm, nrow = length(x_norm), ncol = length(x_norm))
  x2 + t(x2) - 2 * tcrossprod(x_mat)
}

#' @keywords internal
rbf_kernel_matrix <- function(x_mat, scaling) {
  sq_dist <- squared_distance_matrix(x_mat)
  k_mat <- exp(-sq_dist / scaling)
  warn_vanishing_kernel(k_mat)
  k_mat
}

#' @keywords internal
max_offdiag_from_matrix <- function(k_mat) {
  if (!is.matrix(k_mat) || nrow(k_mat) < 2) {
    return(Inf)
  }
  vals <- k_mat[lower.tri(k_mat, diag = FALSE)]
  if (length(vals) == 0) {
    return(Inf)
  }
  max(vals, na.rm = TRUE)
}

#' @keywords internal
compute_v_matrix <- function(X, grads, scaling) {
  n <- nrow(X)
  d <- ncol(X)

  sq_dist <- squared_distance_matrix(X)
  k_mat <- exp(-sq_dist / scaling)
  warn_vanishing_kernel(k_mat)

  pairwise_scores <- tcrossprod(grads)

  score_dot_x <- rowSums(grads * X)
  score_x_cross <- tcrossprod(grads, X)

  term_b <- (2 / scaling) * k_mat *
    (matrix(score_dot_x, n, n) - score_x_cross)
  term_c <- (2 / scaling) * k_mat *
    (matrix(score_dot_x, n, n, byrow = TRUE) - t(score_x_cross))

  term_d <- k_mat * ((2 * d / scaling) - (4 * sq_dist / (scaling^2)))

  pairwise_scores * k_mat + term_b + term_c + term_d
}

#' @keywords internal
warn_vanishing_kernel <- function(k_mat, threshold = 1e-12) {
  if (!is.matrix(k_mat) || length(k_mat) == 0 || nrow(k_mat) < 2) {
    return(invisible(NULL))
  }

  offdiag_max <- max_offdiag_from_matrix(k_mat)
  if (is.finite(offdiag_max) && offdiag_max < threshold) {
    warning(
      "Vanishing off-diagonal kernel detected: bandwidth may be too small and samples become effectively isolated.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @keywords internal
compute_v_matrix_single_dimension <- function(X, grads, scaling, dim_index, k_mat = NULL) {
  if (is.null(k_mat)) {
    x_dim_mat <- matrix(X[, dim_index], ncol = 1)
    k_mat <- rbf_kernel_matrix(x_dim_mat, scaling)
  }

  x_dim <- X[, dim_index]
  g_dim <- grads[, dim_index]

  diffs <- outer(x_dim, x_dim, "-")
  sq_diff <- diffs^2

  term_a <- tcrossprod(g_dim) * k_mat
  g1k <- (-2 / scaling) * k_mat * diffs
  g2k <- -g1k
  term_b <- sweep(g1k, 2, g_dim, "*")
  term_c <- sweep(g2k, 1, g_dim, "*")
  term_d <- 2 * k_mat * (scaling - 2 * sq_diff) / (scaling^2)

  term_a + term_b + term_c + term_d
}

#' @keywords internal
compute_v_block <- function(X_i, X_j, g_i, g_j, scaling) {
  d <- ncol(X_i)
  x_i_norm <- rowSums(X_i * X_i)
  x_j_norm <- rowSums(X_j * X_j)
  sq_dist <- outer(x_i_norm, x_j_norm, "+") - 2 * tcrossprod(X_i, X_j)
  k_mat <- exp(-sq_dist / scaling)

  pairwise_scores <- tcrossprod(g_i, g_j)

  score_dot_i <- rowSums(g_i * X_i)
  score_dot_j <- rowSums(g_j * X_j)
  score_x_cross_ij <- tcrossprod(g_i, X_j)
  score_x_cross_ji <- tcrossprod(g_j, X_i)

  term_b <- (2 / scaling) * k_mat *
    (matrix(score_dot_i, nrow(X_i), nrow(X_j)) - score_x_cross_ij)
  term_c <- (2 / scaling) * k_mat *
    (matrix(score_dot_j, nrow(X_i), nrow(X_j), byrow = TRUE) - t(score_x_cross_ji))

  term_d <- k_mat * ((2 * d / scaling) - (4 * sq_dist / (scaling^2)))

  list(
    u_block = pairwise_scores * k_mat + term_b + term_c + term_d,
    k_block = k_mat
  )
}

#' @keywords internal
compute_v_block_1d <- function(x_i, x_j, g_i, g_j, scaling) {
  diffs <- outer(x_i, x_j, "-")
  sq_diff <- diffs^2
  k_mat <- exp(-sq_diff / scaling)

  term_a <- tcrossprod(g_i, g_j) * k_mat
  g1k <- (-2 / scaling) * k_mat * diffs
  g2k <- -g1k
  term_b <- sweep(g1k, 2, g_j, "*")
  term_c <- sweep(g2k, 1, g_i, "*")
  term_d <- 2 * k_mat * (scaling - 2 * sq_diff) / (scaling^2)

  list(
    u_block = term_a + term_b + term_c + term_d,
    k_block = k_mat
  )
}

#' @keywords internal
generate_bootstrap_weights <- function(n, nboot, boot_method, change_prob = NULL) {
  if (boot_method == "rademacher") {
    matrix(ifelse(stats::rnorm(n * nboot) >= 0, 1, -1), nrow = n, ncol = nboot)
  } else if (boot_method == "markov") {
    replicate(nboot, simulatepm(n, change_prob))
  } else {
    stop("Unknown boot_method")
  }
}

#' @keywords internal
compute_v_stat_bootstrap_blocked <- function(X, grads, scaling, boot_method, nboot, change_prob, block_size, dim_index = NULL) {
  n <- nrow(X)
  weights <- generate_bootstrap_weights(n, nboot, boot_method, change_prob)

  boot_stats <- rep(0, nboot)
  sum_u <- 0
  max_offdiag_kernel <- -Inf

  i_starts <- seq(1, n, by = block_size)
  j_starts <- seq(1, n, by = block_size)

  for (i_start in i_starts) {
    i_end <- min(i_start + block_size - 1, n)
    ii <- i_start:i_end

    if (is.null(dim_index)) {
      X_i <- X[ii, , drop = FALSE]
      g_i <- grads[ii, , drop = FALSE]
    } else {
      x_i <- X[ii, dim_index]
      g_i <- grads[ii, dim_index]
    }
    w_i <- weights[ii, , drop = FALSE]

    for (j_start in j_starts) {
      j_end <- min(j_start + block_size - 1, n)
      jj <- j_start:j_end

      if (is.null(dim_index)) {
        X_j <- X[jj, , drop = FALSE]
        g_j <- grads[jj, , drop = FALSE]
        block <- compute_v_block(X_i, X_j, g_i, g_j, scaling)
      } else {
        x_j <- X[jj, dim_index]
        g_j <- grads[jj, dim_index]
        block <- compute_v_block_1d(x_i, x_j, g_i, g_j, scaling)
      }

      sum_u <- sum_u + sum(block$u_block)

      if (i_start == j_start) {
        if (nrow(block$k_block) > 1) {
          local_max <- max(block$k_block[lower.tri(block$k_block, diag = FALSE)])
          max_offdiag_kernel <- max(max_offdiag_kernel, local_max)
        }
      } else {
        max_offdiag_kernel <- max(max_offdiag_kernel, max(block$k_block))
      }

      w_j <- weights[jj, , drop = FALSE]
      boot_stats <- boot_stats + colSums((block$u_block %*% w_j) * w_i) / n
    }
  }

  if (is.finite(max_offdiag_kernel) && max_offdiag_kernel < 1e-12) {
    warning(
      "Vanishing off-diagonal kernel detected: bandwidth may be too small and samples become effectively isolated.",
      call. = FALSE
    )
  }

  stat_value <- sum_u / n
  list(
    stat_value = stat_value,
    p_value = mean(boot_stats > stat_value),
    boot_stats = as.numeric(boot_stats),
    max_offdiag_kernel = max_offdiag_kernel
  )
}

#' @keywords internal
perform_bootstrap <- function(U_matrix, boot_method, nboot, change_prob = NULL) {
  n <- nrow(U_matrix)
  stat <- n * mean(U_matrix)

  w_mat <- generate_bootstrap_weights(n, nboot, boot_method, change_prob)
  boot_stats <- colSums((U_matrix %*% w_mat) * w_mat) / n

  list(
    p_value = mean(boot_stats > stat),
    boot_stats = as.numeric(boot_stats)
  )
}
