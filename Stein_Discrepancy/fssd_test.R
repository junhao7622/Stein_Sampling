# fssd_test.R

#' Finite Set Stein Discrepancy (FSSD) Linear Time Test
#'
#' Computes a linear-time Stein goodness-of-fit test using random test
#' frequencies and a Gaussian test function family.
#'
#' @param X Numeric vector or matrix of samples. A vector is treated as `n x 1`.
#' @param score_function Function taking `X` and returning score values with shape
#' `n x d` (or length `n` when `d = 1`).
#' @param num_random_frequencies Number of random frequencies used by FSSD.
#' @param scaling Numeric length-2 vector giving the random scale range.
#'
#' @return An object of class `htest` containing the FSSD test result and
#' diagnostic information.
#'
#' @details
#' Let `normal_under_null` be the `n x (d * J)` matrix of Stein features,
#' where `d = ncol(X)` and `J = num_random_frequencies`. For each frequency,
#' a full `n x d` Stein feature matrix is computed and horizontally stacked.
#' The test statistic is the Mahalanobis-type quadratic form of the column means
#' of this matrix, scaled by sample size. A chi-squared reference distribution
#' with `d * J` degrees of freedom is used.
#'
#' @examples
#' model <- gmm(nComp = 3, d = 2)
#' X <- rgmm(model, n = 100)
#' grad_log_prob <- get_score_evaluator(model)
#' out <- fssd_test(X, grad_log_prob, num_random_frequencies = 5)
#' out$p.value
#'
#' @export
fssd_test <- function(X, score_function,
                      num_random_frequencies = 5,
                      scaling = c(1.0, 10.0)) {
  data_name <- deparse(substitute(X))

  x_mat <- as_samples_matrix_fssd(X)
  n <- nrow(x_mat)
  d <- ncol(x_mat)

  if (!is.function(score_function)) {
    stop("score_function must be a function")
  }
  if (!is.numeric(num_random_frequencies) || length(num_random_frequencies) != 1 ||
      num_random_frequencies <= 0) {
    stop("num_random_frequencies must be a positive scalar")
  }
  num_random_frequencies <- as.integer(num_random_frequencies)

  if (!is.numeric(scaling) || length(scaling) != 2 || any(!is.finite(scaling)) ||
      any(scaling <= 0) || scaling[1] > scaling[2]) {
    stop("scaling must be a positive numeric vector of length 2 with scaling[1] <= scaling[2]")
  }

  grads <- as_gradient_matrix_fssd(score_function, x_mat)

  normal_under_null <- matrix(0.0, nrow = n, ncol = num_random_frequencies * d)
  for (freq_index in seq_len(num_random_frequencies)) {
    random_frequency <- stats::rnorm(d)
    random_scale <- stats::runif(1, scaling[1], scaling[2])

    start_col <- (freq_index - 1) * d + 1
    end_col <- freq_index * d
    normal_under_null[, start_col:end_col] <-
      compute_stein_stat_fssd(x_mat, grads, random_frequency, random_scale)
  }

  maha <- compute_mahalanobis_pvalue(normal_under_null)

  result <- list(
    statistic = c(fssd = maha$statistic),
    p.value = maha$p_value,
    method = "Finite Set Stein Discrepancy (Linear Time)",
    data.name = data_name,
    info = list(
      scaling_range = scaling,
      num_random_frequencies = num_random_frequencies,
      normal_under_null = normal_under_null
    )
  )

  class(result) <- "htest"
  result
}

#' Coerce input samples to an `n x d` numeric matrix
#'
#' @param x Input samples as vector or matrix.
#'
#' @return Numeric matrix with at least two rows.
#' @keywords internal
as_samples_matrix_fssd <- function(x) {
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

#' Validate and coerce score output to an `n x d` matrix
#'
#' @param score_function Score function supplied by user.
#' @param x_mat Sample matrix (`n x d`).
#'
#' @return Numeric matrix of score values matching `x_mat` dimensions.
#' @keywords internal
as_gradient_matrix_fssd <- function(score_function, x_mat) {
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

#' Compute chi-squared p-value from Stein feature matrix
#'
#' @param difference_matrix Matrix of Stein features (`n x J`).
#'
#' @return List with numeric fields `statistic` and `p_value`.
#' @keywords internal
compute_mahalanobis_pvalue <- function(difference_matrix) {
  x <- as.matrix(difference_matrix)
  n <- nrow(x)
  p <- ncol(x)

  mu <- colMeans(x)

  if (p == 1) {
    sigma <- stats::var(x[, 1])
    if (!is.finite(sigma) || sigma <= 0) {
      return(list(statistic = 0.0, p_value = 1.0))
    }
    stat <- as.numeric(n * (mu^2) / sigma)
    p_val <- stats::pchisq(stat, df = 1, lower.tail = FALSE)
    return(list(statistic = stat, p_value = p_val))
  }

  sigma <- stats::cov(x)
  sigma_inv <- tryCatch(
    solve(sigma),
    error = function(e) {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("MASS is required for pseudo-inverse fallback (MASS::ginv)")
      }
      MASS::ginv(sigma)
    }
  )

  stat <- as.numeric(n * (t(mu) %*% sigma_inv %*% mu))
  if (!is.finite(stat) || stat < 0) {
    return(list(statistic = 0.0, p_value = 1.0))
  }

  p_val <- stats::pchisq(stat, df = p, lower.tail = FALSE)
  list(statistic = stat, p_value = p_val)
}

#' Compute one FSSD Stein feature column (length `n`)
#'
#' @param X Sample matrix (`n x d`).
#' @param grads Score matrix (`n x d`).
#' @param random_frequency Numeric vector of length `d`.
#' @param random_scale Positive scalar bandwidth-like scale.
#'
#' @return Numeric matrix of shape `n x d` containing Stein feature values.
#' @keywords internal
compute_stein_stat_fssd <- function(X, grads, random_frequency, random_scale) {
  n <- nrow(X)
  d <- ncol(X)

  if (length(random_frequency) != d) {
    stop("random_frequency must have length equal to ncol(X)")
  }
  if (!is.numeric(random_scale) || length(random_scale) != 1 || random_scale <= 0) {
    stop("random_scale must be a positive scalar")
  }

  delta <- X - matrix(random_frequency, nrow = n, ncol = d, byrow = TRUE)
  sq_norm <- rowSums(delta * delta)
  mean_embedding <- exp(-sq_norm / random_scale)

  test_function_grad <- -(2 / random_scale) * delta * mean_embedding
  grads * mean_embedding + test_function_grad
}
