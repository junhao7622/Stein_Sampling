# ksd_u_test.R

if(file.exists("stein_helpers.R")){
  source("stein_helpers.R")
}

if (!exists("as_samples_matrix", mode = "function")) {
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
}

if (!exists("as_gradient_matrix", mode = "function")) {
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

    if (!is.numeric(grads) || nrow(grads) != n || ncol(grads) != d) {
      stop("score_function output shape must match X (n x d)")
    }

    grads
  }
}

#' Kernel Stein Discrepancy GOF Test (U-statistics)
#'
#' Performs a goodness-of-fit test using a Stein discrepancy with U-statistics
#' and bootstrap calibration.
#'
#' @param X Numeric vector or matrix of samples. A vector is treated as `n x 1`.
#' @param score_function Function of `X` returning score values with shape `n x d`
#'   (or a length-`n` vector when `d = 1`).
#' @param boot_method Bootstrap method for degenerate U-statistics.
#'   Only `"multinomial_centered"` is supported.
#' @param scaling Positive kernel scaling parameter. If `NULL`, uses the median
#'   pairwise squared-distance median heuristic.
#' @param nboot Number of bootstrap samples.
#' @param median_max_samples Maximum sample size used by median heuristic.
#' @param median_use_sampling If `TRUE`, median heuristic uses subsampling above
#'   `median_max_samples`.
#' @param median_seed Random seed used by heuristic sampling.
#' @param block_size Optional block size for matrix-free computation.
#' @param block_threshold Auto-enable block mode when `n > block_threshold`.
#' @param kernel Kernel choice: `"rbf"`, `"imq"`, or custom list with fields
#'   `name`, `eval(r2, scaling)`, `grad_factor(r2, scaling)`,
#'   `trace_mixed(r2, scaling, d)`.
#' @param imq_beta Exponent for IMQ kernel `(scaling + ||x-y||^2)^imq_beta`.
#'   Must be negative.
#'
#' @return An object of class `htest` with test statistic, p-value,
#' bootstrap samples, and diagnostics in `info`.
#'
#' @export
ksd_u_test <- function(X, score_function,
                       boot_method = "multinomial_centered",
                       scaling = NULL,
                       nboot = 1000,
                       median_max_samples = 2000,
                       median_use_sampling = TRUE,
                       median_seed = 123,
                       block_size = NULL,
                       block_threshold = 5000,
                       kernel = c("rbf", "imq"),
                       imq_beta = -0.5) {
  data_name <- deparse(substitute(X))
  if (!is.character(boot_method) || length(boot_method) != 1) {
    stop("boot_method must be 'multinomial_centered'")
  }
  if (identical(boot_method, "weighted")) {
    warning("'weighted' is deprecated; use 'multinomial_centered'.", call. = FALSE)
    boot_method <- "multinomial_centered"
  }
  if (!identical(boot_method, "multinomial_centered")) {
    stop("Only boot_method = 'multinomial_centered' is supported for degenerate U-statistics")
  }

  x_mat <- as_samples_matrix(X)
  n <- nrow(x_mat)
  d <- ncol(x_mat)

  if (!is.function(score_function)) {
    stop("score_function must be a function")
  }

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

  kernel_spec <- build_stein_kernel(kernel, scaling = scaling, imq_beta = imq_beta)
  block_info <- resolve_u_block_settings(n, block_size, block_threshold)

  if (block_info$use_block_mode) {
    blocked <- compute_u_stat_bootstrap_blocked(
      x_mat = x_mat,
      grads = grads,
      kernel_spec = kernel_spec,
      scaling = scaling,
      d = d,
      nboot = nboot,
      block_size = block_info$block_size
    )

    ksd_u <- blocked$stat_value
    boot_stats <- blocked$boot_stats
    p_value <- mean(boot_stats >= ksd_u)

    result <- list(
      statistic = c(ksd_u = ksd_u),
      p.value = p_value,
      method = paste("Kernel Stein Discrepancy (U-statistics)", "-", kernel_spec$name),
      data.name = data_name,
      bootstrap_samples = boot_stats,
      info = list(
        scaling = scaling,
        kernel = kernel_spec$name,
        imq_beta = if (kernel_spec$name == "imq") imq_beta else NA_real_,
        M_matrix = NULL,
        nboot = nboot,
        boot_method = boot_method,
        block_mode = TRUE,
        block_size = block_info$block_size,
        max_offdiag_kernel = blocked$max_offdiag_kernel
      )
    )
  } else {
    m_matrix <- compute_u_matrix(x_mat, grads, scaling, d, kernel_spec)
    m_u_matrix <- m_matrix - diag(diag(m_matrix))
    ksd_u <- sum(m_u_matrix) / (n * (n - 1))

    boot_stats <- perform_u_bootstrap(m_u_matrix, nboot)
    p_value <- mean(boot_stats >= ksd_u)

    result <- list(
      statistic = c(ksd_u = ksd_u),
      p.value = p_value,
      method = paste("Kernel Stein Discrepancy (U-statistics)", "-", kernel_spec$name),
      data.name = data_name,
      bootstrap_samples = boot_stats,
      info = list(
        scaling = scaling,
        kernel = kernel_spec$name,
        imq_beta = if (kernel_spec$name == "imq") imq_beta else NA_real_,
        M_matrix = m_matrix,
        nboot = nboot,
        boot_method = boot_method,
        block_mode = FALSE,
        block_size = NA_integer_,
        max_offdiag_kernel = max_offdiag_from_kernel_matrix(kernel_spec$eval(squared_distance_matrix_u(x_mat), scaling))
      )
    )
  }

  class(result) <- "htest"
  result
}

#' @keywords internal
resolve_u_block_settings <- function(n, block_size, block_threshold) {
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
build_stein_kernel <- function(kernel, scaling, imq_beta = -0.5) {
  if (is.list(kernel)) {
    required <- c("name", "eval", "grad_factor", "trace_mixed")
    if (!all(required %in% names(kernel))) {
      stop("Custom kernel list must contain: name, eval, grad_factor, trace_mixed")
    }
    if (!is.function(kernel$eval) || !is.function(kernel$grad_factor) || !is.function(kernel$trace_mixed)) {
      stop("Custom kernel entries eval, grad_factor, trace_mixed must be functions")
    }
    return(kernel)
  }

  kernel_name <- match.arg(kernel, c("rbf", "imq"))

  if (kernel_name == "rbf") {
    return(list(
      name = "rbf",
      eval = function(r2, scaling) exp(-r2 / scaling),
      grad_factor = function(r2, scaling) (-2 / scaling) * exp(-r2 / scaling),
      trace_mixed = function(r2, scaling, d) {
        k <- exp(-r2 / scaling)
        k * ((2 * d / scaling) - (4 * r2 / (scaling^2)))
      }
    ))
  }

  if (!is.numeric(imq_beta) || length(imq_beta) != 1 || !is.finite(imq_beta) || imq_beta >= 0) {
    stop("imq_beta must be a finite negative scalar")
  }

  list(
    name = "imq",
    eval = function(r2, scaling) (scaling + r2)^imq_beta,
    grad_factor = function(r2, scaling) 2 * imq_beta * (scaling + r2)^(imq_beta - 1),
    trace_mixed = function(r2, scaling, d) {
      s <- scaling + r2
      -4 * imq_beta * (imq_beta - 1) * r2 * s^(imq_beta - 2) -
        2 * imq_beta * d * s^(imq_beta - 1)
    }
  )
}

#' @keywords internal
squared_distance_matrix_u <- function(x_mat) {
  x_norm <- rowSums(x_mat * x_mat)
  x2 <- matrix(x_norm, nrow = length(x_norm), ncol = length(x_norm))
  x2 + t(x2) - 2 * tcrossprod(x_mat)
}

#' @keywords internal
max_offdiag_from_kernel_matrix <- function(k_mat) {
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
warn_vanishing_kernel <- function(k_mat, threshold = 1e-12) {
  if (!is.matrix(k_mat) || length(k_mat) == 0 || nrow(k_mat) < 2) {
    return(invisible(NULL))
  }

  offdiag_max <- max_offdiag_from_kernel_matrix(k_mat)
  if (is.finite(offdiag_max) && offdiag_max < threshold) {
    warning(
      "Vanishing off-diagonal kernel detected: bandwidth may be too small and samples become effectively isolated.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @keywords internal
compute_u_matrix <- function(x_mat, grads, scaling, d, kernel_spec) {
  n <- nrow(x_mat)
  pairwise_sq_dist <- squared_distance_matrix_u(x_mat)

  k_xy <- kernel_spec$eval(pairwise_sq_dist, scaling)
  warn_vanishing_kernel(k_xy)

  grad_factor <- kernel_spec$grad_factor(pairwise_sq_dist, scaling)
  score_dot <- rowSums(x_mat * grads)
  score_x_cross <- tcrossprod(grads, x_mat)

  score_x_dy <- -grad_factor * (matrix(score_dot, n, n) - score_x_cross)
  dx_score_y <- grad_factor * (t(score_x_cross) - matrix(score_dot, n, n, byrow = TRUE))

  trace_xy <- kernel_spec$trace_mixed(pairwise_sq_dist, scaling, d)

  tcrossprod(grads) * k_xy + score_x_dy + dx_score_y + trace_xy
}

#' @keywords internal
compute_u_block <- function(X_i, X_j, g_i, g_j, scaling, d, kernel_spec) {
  ni <- nrow(X_i)
  nj <- nrow(X_j)

  x_i_norm <- rowSums(X_i * X_i)
  x_j_norm <- rowSums(X_j * X_j)
  sq_dist <- outer(x_i_norm, x_j_norm, "+") - 2 * tcrossprod(X_i, X_j)

  k_block <- kernel_spec$eval(sq_dist, scaling)
  grad_factor <- kernel_spec$grad_factor(sq_dist, scaling)

  score_dot_i <- rowSums(X_i * g_i)
  score_dot_j <- rowSums(X_j * g_j)
  score_x_cross_ij <- tcrossprod(g_i, X_j)
  score_x_cross_ji <- tcrossprod(g_j, X_i)

  term_a <- tcrossprod(g_i, g_j) * k_block
  term_b <- -grad_factor * (matrix(score_dot_i, ni, nj) - score_x_cross_ij)
  term_c <- grad_factor * (t(score_x_cross_ji) - matrix(score_dot_j, ni, nj, byrow = TRUE))
  term_d <- kernel_spec$trace_mixed(sq_dist, scaling, d)

  list(
    m_block = term_a + term_b + term_c + term_d,
    k_block = k_block
  )
}

#' @keywords internal
compute_u_stat_bootstrap_blocked <- function(x_mat, grads, kernel_spec, scaling, d, nboot, block_size) {
  n <- nrow(x_mat)
  weight_mat <- stats::rmultinom(nboot, size = n, prob = rep(1 / n, n)) / n
  vec_mat <- weight_mat - (1 / n)

  i_starts <- seq(1, n, by = block_size)
  j_starts <- seq(1, n, by = block_size)

  sum_offdiag <- 0
  quad_sum <- rep(0, nboot)
  max_offdiag_kernel <- -Inf

  for (i_start in i_starts) {
    i_end <- min(i_start + block_size - 1, n)
    ii <- i_start:i_end

    X_i <- x_mat[ii, , drop = FALSE]
    g_i <- grads[ii, , drop = FALSE]
    v_i <- vec_mat[ii, , drop = FALSE]

    for (j_start in j_starts) {
      j_end <- min(j_start + block_size - 1, n)
      jj <- j_start:j_end

      X_j <- x_mat[jj, , drop = FALSE]
      g_j <- grads[jj, , drop = FALSE]
      v_j <- vec_mat[jj, , drop = FALSE]

      block <- compute_u_block(X_i, X_j, g_i, g_j, scaling, d, kernel_spec)
      m_block <- block$m_block

      if (i_start == j_start) {
        diag(m_block) <- 0
        if (nrow(block$k_block) > 1) {
          local_max <- max(block$k_block[lower.tri(block$k_block, diag = FALSE)])
          max_offdiag_kernel <- max(max_offdiag_kernel, local_max)
        }
      } else {
        max_offdiag_kernel <- max(max_offdiag_kernel, max(block$k_block))
      }

      sum_offdiag <- sum_offdiag + sum(m_block)
      quad_sum <- quad_sum + colSums((m_block %*% v_j) * v_i)
    }
  }

  if (is.finite(max_offdiag_kernel) && max_offdiag_kernel < 1e-12) {
    warning(
      "Vanishing off-diagonal kernel detected: bandwidth may be too small and samples become effectively isolated.",
      call. = FALSE
    )
  }

  stat_value <- sum_offdiag / (n * (n - 1))
  boot_stats <- quad_sum

  list(
    stat_value = stat_value,
    boot_stats = as.numeric(boot_stats),
    max_offdiag_kernel = max_offdiag_kernel
  )
}

#' @keywords internal
perform_u_bootstrap <- function(m_u_matrix, nboot) {
  n <- nrow(m_u_matrix)
  weight_mat <- stats::rmultinom(nboot, size = n, prob = rep(1 / n, n)) / n
  centered_mat <- weight_mat - (1 / n)
  boot_stats <- colSums((m_u_matrix %*% centered_mat) * centered_mat)

  as.numeric(boot_stats)
}

# Backward-compatible aliases
ksd_u_stats <- ksd_u_test
KSD <- ksd_u_test
