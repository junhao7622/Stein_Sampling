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
#' Performs a goodness-of-fit test using a Gaussian-kernel Stein discrepancy
#' with U-statistics and bootstrap calibration.
#'
#' @param X Numeric vector or matrix of samples. A vector is treated as `n x 1`.
#' @param score_function Function of `X` returning score values with shape `n x d`
#'   (or a length-`n` vector when `d = 1`).
#' @param boot_method One of `"weighted"` or `"efron"`.
#' @param scaling Positive kernel scaling parameter. If `NULL`, uses the median
#'   pairwise-distance heuristic and converts it to Gaussian-kernel scale.
#' @param nboot Number of bootstrap samples.
#'
#' @return An object of class `htest` with test statistic, p-value,
#' bootstrap samples, and diagnostics in `info`.
#'
#' The returned `info` list contains:
#' \itemize{
#' \item `bandwidth`: Gaussian kernel bandwidth used in the test.
#' \item `M_matrix`: Full Stein matrix before diagonal removal.
#' \item `nboot`: Number of bootstrap replications.
#' \item `boot_method`: Bootstrap strategy used for calibration.
#' }
#'
#' @details
#' Let \eqn{M} denote the pairwise Stein kernel matrix and \eqn{M_u} be its
#' off-diagonal counterpart (diagonal set to zero). The reported U-statistic is
#' \eqn{\widehat{\mathrm{KSD}}_u = \sum_{i \neq j} M_{ij} / (n(n-1))}. P-values
#' are computed by bootstrap from \eqn{M_u} using either weighted multinomial
#' resampling or Efron index-resampling.
#'
#' @examples
#' model <- gmm(nComp = 3, d = 2)
#' X <- rgmm(model, n = 80)
#' grad_log_prob <- get_score_evaluator(model)
#'
#' out <- ksd_u_test(X, grad_log_prob, boot_method = "weighted", nboot = 300)
#' out$statistic
#' out$p.value
#'
#' @export
ksd_u_test <- function(X, score_function,
                       boot_method = c("weighted", "efron"),
                       scaling = NULL,
                       nboot = 1000) {
  data_name <- deparse(substitute(X))
  boot_method <- match.arg(boot_method)

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

  grads <- as_gradient_matrix(score_function, x_mat)

  if (is.null(scaling)) {
    scaling <- sqrt(0.5 * find_median_distance(x_mat))
  }
  if (!is.numeric(scaling) || length(scaling) != 1 || !is.finite(scaling) || scaling <= 0) {
    stop("scaling must be a positive scalar")
  }

  m_matrix <- compute_u_matrix(x_mat, grads, scaling, d)
  m_u_matrix <- m_matrix - diag(diag(m_matrix))
  ksd_u <- sum(m_u_matrix) / (n * (n - 1))

  boot_stats <- perform_u_bootstrap(m_u_matrix, ksd_u, boot_method, nboot)
  p_value <- mean(boot_stats >= ksd_u)

  result <- list(
    statistic = c(ksd_u = ksd_u),
    p.value = p_value,
    method = "Kernel Stein Discrepancy (U-statistics)",
    data.name = data_name,
    bootstrap_samples = boot_stats,
    info = list(
      bandwidth = scaling,
      M_matrix = m_matrix,
      nboot = nboot,
      boot_method = boot_method
    )
  )
  class(result) <- "htest"
  result
}

#' @keywords internal
compute_u_matrix <- function(x_mat, grads, bandwidth, d) {
  xy <- x_mat %*% t(x_mat)
  x_sq <- matrix(rowSums(x_mat * x_mat), ncol = 1)
  x_sq_expand <- rep_mat(x_sq, 1, nrow(x_mat))

  pairwise_sq_dist <- x_sq_expand + t(x_sq_expand) - 2 * xy
  k_xy <- exp(-pairwise_sq_dist / (2 * bandwidth^2))

  score_dot <- rowSums(x_mat * grads)
  score_x_dy <- -(grads %*% t(x_mat) - rep_mat(score_dot, 1, nrow(x_mat))) / (bandwidth^2)
  dx_score_y <- t(score_x_dy)
  dx_dy <- -pairwise_sq_dist / (bandwidth^4) + d / (bandwidth^2)

  (grads %*% t(grads) + score_x_dy + dx_score_y + dx_dy) * k_xy
}

#' @keywords internal
perform_u_bootstrap <- function(m_u_matrix, stat_value, boot_method, nboot) {
  n <- nrow(m_u_matrix)
  boot_stats <- numeric(nboot)

  if (boot_method == "weighted") {
    for (boot_index in seq_len(nboot)) {
      boot_weights <- as.numeric(stats::rmultinom(1, size = n, prob = rep(1 / n, n))) / n
      centered_weights <- boot_weights - (1 / n)
      boot_stats[boot_index] <- as.numeric(t(centered_weights) %*% m_u_matrix %*% centered_weights)
    }
  } else if (boot_method == "efron") {
    for (boot_index in seq_len(nboot)) {
      boot_index_sample <- sample.int(n, size = n, replace = TRUE)
      boot_stats[boot_index] <- sum(m_u_matrix[boot_index_sample, boot_index_sample]) / (n * (n - 1))
    }
  } else {
    stop("Unknown boot_method")
  }

  boot_stats
}

# Backward-compatible aliases
ksd_u_stats <- ksd_u_test
KSD <- ksd_u_test
