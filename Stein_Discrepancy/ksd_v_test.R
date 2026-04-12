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
#' @param test_type One of `"aggregated"` or `"per_dimension"`.
#' @param boot_method One of `"rademacher"` or `"markov"`.
#' @param scaling Positive kernel scaling parameter. If `NULL`, uses the median
#'   squared-distance heuristic.
#' @param nboot Number of bootstrap samples.
#' @param change_prob Markov sign-flip probability when `boot_method = "markov"`.
#'
#' @return An object of class `htest` with test statistic(s), p-value information,
#' bootstrap samples, and diagnostic matrices in `info`.
#'
#' @export
ksd_v_test <- function(X, score_function,
                       test_type = c("aggregated", "per_dimension"),
                       boot_method = c("rademacher", "markov"),
                       scaling = NULL,
                       nboot = 1000,
                       change_prob = 0.1) {
  data_name <- deparse(substitute(X))
  test_type <- match.arg(test_type)
  boot_method <- match.arg(boot_method)

  x_mat <- as_samples_matrix(X)
  n <- nrow(x_mat)
  d <- ncol(x_mat)

  if (!is.numeric(nboot) || length(nboot) != 1 || nboot <= 0) {
    stop("nboot must be a positive scalar")
  }
  nboot <- as.integer(nboot)

  if (!is.numeric(change_prob) || length(change_prob) != 1 ||
      change_prob < 0 || change_prob > 1) {
    stop("change_prob must be in [0, 1]")
  }

  grads <- as_gradient_matrix(score_function, x_mat)

  if (is.null(scaling)) {
    scaling <- find_median_distance(x_mat)
  }

  if (!is.numeric(scaling) || length(scaling) != 1 || scaling <= 0) {
    stop("scaling must be a positive scalar")
  }

  if (test_type == "aggregated") {
    u_mat <- compute_v_matrix(x_mat, grads, scaling)
    stat_value <- n * mean(u_mat)
    boot <- perform_bootstrap(u_mat, boot_method, nboot, change_prob)

    result <- list(
      statistic = c(ksd_v = stat_value),
      p.value = boot$p_value,
      method = paste("Kernel Stein Discrepancy (V-statistics) -", test_type),
      data.name = data_name,
      bootstrap_samples = boot$boot_stats,
      info = list(
        bandwidth = scaling,
        U_matrix = u_mat,
        nboot = nboot,
        boot_method = boot_method
      )
    )
  } else {
    k_mat <- rbf_kernel_matrix(x_mat, scaling)
    stats <- numeric(d)
    p_vals <- numeric(d)
    boot_mat <- matrix(NA_real_, nrow = nboot, ncol = d)
    u_mats <- vector("list", d)

    for (dim_index in seq_len(d)) {
      u_dim <- compute_v_matrix_single_dimension(x_mat, grads, scaling, dim_index, k_mat)
      stats[dim_index] <- n * mean(u_dim)
      boot <- perform_bootstrap(u_dim, boot_method, nboot, change_prob)
      p_vals[dim_index] <- boot$p_value
      boot_mat[, dim_index] <- boot$boot_stats
      u_mats[[dim_index]] <- u_dim
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
        bandwidth = scaling,
        U_matrices = u_mats,
        raw_p_values = p_vals,
        nboot = nboot,
        boot_method = boot_method
      )
    )
  }

  class(result) <- "htest"
  result
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
  x2 + t(x2) - 2 * (x_mat %*% t(x_mat))
}

#' @keywords internal
rbf_kernel_matrix <- function(x_mat, scaling) {
  sq_dist <- squared_distance_matrix(x_mat)
  exp(-sq_dist / scaling)
}

#' @keywords internal
compute_v_matrix <- function(X, grads, scaling) {
  n <- nrow(X)
  d <- ncol(X)

  sq_dist <- squared_distance_matrix(X)
  k_mat <- exp(-sq_dist / scaling)

  pairwise_scores <- grads %*% t(grads)

  score_dot_x <- rowSums(grads * X)
  score_x_cross <- grads %*% t(X)

  term_b <- (2 / scaling) * k_mat *
    (matrix(score_dot_x, n, n) - score_x_cross)
  term_c <- (2 / scaling) * k_mat *
    (matrix(score_dot_x, n, n, byrow = TRUE) - t(score_x_cross))

  term_d <- k_mat * ((2 * d / scaling) - (4 * sq_dist / (scaling^2)))

  pairwise_scores * k_mat + term_b + term_c + term_d
}

#' @keywords internal
compute_v_matrix_single_dimension <- function(X, grads, scaling, dim_index, k_mat = NULL) {
  if (is.null(k_mat)) {
    k_mat <- rbf_kernel_matrix(X, scaling)
  }

  x_dim <- X[, dim_index]
  g_dim <- grads[, dim_index]

  diffs <- outer(x_dim, x_dim, "-")
  sq_diff <- diffs^2

  term_a <- outer(g_dim, g_dim) * k_mat
  g1k <- (-2 / scaling) * k_mat * diffs
  g2k <- -g1k
  term_b <- sweep(g1k, 2, g_dim, "*")
  term_c <- sweep(g2k, 1, g_dim, "*")
  term_d <- 2 * k_mat * (scaling - 2 * sq_diff) / (scaling^2)

  term_a + term_b + term_c + term_d
}

#' @keywords internal
perform_bootstrap <- function(U_matrix, boot_method, nboot, change_prob = NULL) {
  n <- nrow(U_matrix)
  stat <- n * mean(U_matrix)

  draw_weights <- switch(
    boot_method,
    rademacher = function() ifelse(rnorm(n) >= 0, 1, -1),
    markov = function() simulatepm(n, change_prob),
    stop("Unknown boot_method")
  )

  boot_stats <- replicate(nboot, {
    w <- draw_weights()
    n * mean(U_matrix * tcrossprod(w))
  })

  list(
    p_value = mean(boot_stats > stat),
    boot_stats = as.numeric(boot_stats)
  )
}
