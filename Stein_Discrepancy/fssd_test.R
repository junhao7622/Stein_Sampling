# fssd_test.R

#' FSSD htclass entry point
#'
#' @param X Numeric vector or matrix.
#' @param score_function Score function.
#' @param scaling Optional legacy parameter kept for alignment with KSD tests.
#'   It is ignored in optimized FSSD.
#' @param nboot Number of null simulations (alias of `n_simulations`).
#' @param num_random_frequencies Number of test locations (mapped to `J`).
#' @param n_simulations Deprecated alias of `nboot`.
#' @param ... Additional arguments to `fssd_test_optimized`.
#'
#' @return An object of class `htest`.
#' @export
fssd_test <- function(X,
                      score_function,
                      scaling = NULL,
                      nboot = 2000,
                      num_random_frequencies = 5,
                      n_simulations = NULL,
                      ...) {
  if (!is.null(n_simulations)) {
    if (!missing(nboot) && !identical(as.numeric(n_simulations), as.numeric(nboot))) {
      warning("Both 'nboot' and deprecated 'n_simulations' were provided; using 'n_simulations'.", call. = FALSE)
    }
    nboot <- n_simulations
  }

  fssd_test_optimized(
    X = X,
    score_function = score_function,
    scaling = scaling,
    nboot = nboot,
    J = as.integer(num_random_frequencies),
    ...
  )
}

#' Finite Set Stein Discrepancy (FSSD) Goodness-of-Fit Test
#'
#' Exact FSSD pipeline with data splitting, parameter optimization,
#' unbiased U-statistic evaluation on a held-out test split, and null
#' simulation from covariance eigenvalues.
#'
#' @param X Numeric vector or matrix of samples. A vector is treated as `n x 1`.
#' @param score_function Function taking `X` and returning scores with shape
#' `n x d` (or length `n` when `d = 1`).
#' @param scaling Optional legacy parameter kept for API alignment with KSD
#'   tests. It is ignored in optimized FSSD.
#' @param nboot Number of null simulations (alias of `n_simulations`).
#' @param J Number of test locations.
#' @param train_ratio Fraction of samples used for parameter optimization.
#' @param n_simulations Deprecated alias of `nboot`.
#' @param gamma Stability constant in objective denominator.
#' @param maxit Maximum optimizer iterations.
#' @param seed Optional random seed.
#'
#' @return An object of class `htest`.
#' @export
fssd_test_optimized <- function(X, score_function,
                                scaling = NULL,
                                nboot = 2000,
                                J,
                                train_ratio = 0.5,
                                n_simulations = NULL,
                                gamma = 1e-4,
                                maxit = 200,
                                seed = NULL,
                                ...) {
  data_name <- deparse(substitute(X))

  if (!is.null(scaling)) {
    warning("Argument 'scaling' is ignored in optimized FSSD.", call. = FALSE)
  }

  if (!is.null(n_simulations)) {
    if (!missing(nboot) && !identical(as.numeric(n_simulations), as.numeric(nboot))) {
      warning("Both 'nboot' and deprecated 'n_simulations' were provided; using 'n_simulations'.", call. = FALSE)
    }
    nboot <- n_simulations
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  x_mat <- as_samples_matrix_fssd(X)
  n <- nrow(x_mat)

  if (!is.function(score_function)) {
    stop("score_function must be a function")
  }
  if (!is.numeric(J) || length(J) != 1 || !is.finite(J) || J <= 0) {
    stop("J must be a positive scalar")
  }
  J <- as.integer(J)

  if (!is.numeric(train_ratio) || length(train_ratio) != 1 || !is.finite(train_ratio) ||
      train_ratio <= 0 || train_ratio >= 1) {
    stop("train_ratio must be a scalar in (0, 1)")
  }
  if (!is.numeric(nboot) || length(nboot) != 1 ||
      !is.finite(nboot) || nboot <= 0) {
    stop("nboot must be a positive scalar")
  }
  nboot <- as.integer(nboot)

  if (!is.numeric(gamma) || length(gamma) != 1 || !is.finite(gamma) || gamma <= 0) {
    stop("gamma must be a positive scalar")
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || !is.finite(maxit) || maxit <= 0) {
    stop("maxit must be a positive scalar")
  }
  maxit <- as.integer(maxit)

  n_train <- floor(train_ratio * n)
  n_train <- max(2L, min(n_train, n - 2L))
  idx_train <- sample.int(n, size = n_train, replace = FALSE)
  idx_test <- setdiff(seq_len(n), idx_train)

  x_train <- x_mat[idx_train, , drop = FALSE]
  x_test <- x_mat[idx_test, , drop = FALSE]

  grads_train <- as_gradient_matrix_fssd(score_function, x_train)
  grads_test <- as_gradient_matrix_fssd(score_function, x_test)

  opt <- optimize_fssd(
    X_train = x_train,
    grads_train = grads_train,
    J = J,
    gamma = gamma,
    maxit = maxit
  )

  stat_out <- compute_fssd_stat(
    X_test = x_test,
    grads_test = grads_test,
    V_opt = opt$V_opt,
    sigma2_opt = opt$sigma2_opt
  )

  p_out <- fssd_pvalue(
    tau_matrix = stat_out$tau_matrix,
    fssd_stat = stat_out$fssd_stat,
    n_simulations = nboot
  )

  result <- list(
    statistic = c(fssd = stat_out$fssd_stat),
    p.value = p_out$p_value,
    method = "Finite Set Stein Discrepancy",
    data.name = data_name,
    info = list(
      train_ratio = train_ratio,
      n_train = n_train,
      n_test = nrow(x_test),
      J = J,
      gamma = gamma,
      n_simulations = nboot,
      sigma2_opt = opt$sigma2_opt,
      V_opt = opt$V_opt,
      objective_opt = opt$objective_opt,
      convergence = opt$convergence,
      test_score = p_out$test_score,
      eigenvalues = p_out$eigenvalues
    )
  )
  class(result) <- "htest"
  result
}

#' Optimize FSSD parameters on train split
#'
#' @param X_train Train samples (`n_tr x d`).
#' @param grads_train Train scores (`n_tr x d`).
#' @param J Number of test locations.
#' @param gamma Stability constant.
#' @param maxit Maximum optimizer iterations.
#'
#' @return List with optimized `V_opt` and `sigma2_opt`.
#' @export
optimize_fssd <- function(X_train,
                          grads_train,
                          J,
                          gamma = 1e-4,
                          maxit = 200) {
  x_train <- as_samples_matrix_fssd(X_train)
  grads <- as.matrix(grads_train)
  if (!is.numeric(grads) || nrow(grads) != nrow(x_train) || ncol(grads) != ncol(x_train)) {
    stop("grads_train must be a numeric matrix with the same shape as X_train")
  }

  if (!is.numeric(J) || length(J) != 1 || !is.finite(J) || J <= 0) {
    stop("J must be a positive scalar")
  }
  J <- as.integer(J)

  if (!is.numeric(gamma) || length(gamma) != 1 || !is.finite(gamma) || gamma <= 0) {
    stop("gamma must be a positive scalar")
  }
  if (!is.numeric(maxit) || length(maxit) != 1 || !is.finite(maxit) || maxit <= 0) {
    stop("maxit must be a positive scalar")
  }
  maxit <- as.integer(maxit)

  n <- nrow(x_train)
  d <- ncol(x_train)
  n_params <- d * J + 1L
  if (n < n_params) {
    warning(
      "Warning: The number of optimization parameters (d*J+1) exceeds the training sample size. FSSD-opt is ill-posed and highly prone to overfitting. Consider reducing J or increasing train_ratio.",
      call. = FALSE
    )
  }

  mu <- colMeans(x_train)
  sigma <- safe_covariance(x_train)
  v_init <- sample_gaussian_rows(J, mu, sigma)

  med_dist <- median_pairwise_distance(x_train)
  sigma2_init <- max(med_dist^2, 1e-8)

  par_init <- c(as.numeric(t(v_init)), log(sigma2_init))

  objective_and_grad_from_par <- function(par) {
    v_vec <- par[seq_len(J * d)]
    log_sigma2 <- par[J * d + 1L]

    v_mat <- matrix(v_vec, nrow = J, ncol = d, byrow = TRUE)
    sigma2_k <- exp(log_sigma2)

    tau <- compute_tau(x_train, grads, v_mat, sigma2_k)

    n_local <- nrow(tau)
    p_local <- ncol(tau)
    c_u <- n_local * (n_local - 1)

    sum_tau <- colSums(tau)
    fssd2 <- unbiased_fssd_stat(tau)

    # Use the unbiased U-statistic directly in optimization to avoid
    # dead zones caused by hard clipping at 0.
    grad_f <- (2 / c_u) * (matrix(sum_tau, nrow = n_local, ncol = p_local, byrow = TRUE) - tau)

    mu_tau <- colMeans(tau)
    centered_tau <- sweep(tau, 2, mu_tau, "-")
    sigma_tau <- crossprod(centered_tau) / (n_local - 1)
    sigma_h1_sq_raw <- as.numeric((4 / n_local) * t(mu_tau) %*% sigma_tau %*% mu_tau)
    sigma_h1_sq <- max(sigma_h1_sq_raw, 0)
    sigma_h1 <- sqrt(sigma_h1_sq)
    denom <- sigma_h1 + gamma

    # d(sigma_h1^2)/d tau = 8/n * 1 * (Sigma mu)^T + 8/(n-1) * (C mu) * mu^T
    sigma_mu <- as.numeric(sigma_tau %*% mu_tau)
    c_mu <- as.numeric(centered_tau %*% mu_tau)
    grad_sigma_sq <-
      (8 / n_local) * matrix(sigma_mu, nrow = n_local, ncol = p_local, byrow = TRUE) +
      (8 / (n_local - 1)) * (c_mu %o% mu_tau)
    grad_sigma_sq <- grad_sigma_sq / n_local

    if (sigma_h1 > 1e-12) {
      grad_sigma <- grad_sigma_sq / (2 * sigma_h1)
    } else {
      grad_sigma <- matrix(0, nrow = n_local, ncol = p_local)
    }

    grad_obj_tau <- grad_f / denom - (fssd2 / (denom^2)) * grad_sigma

    scale_fac <- 1 / sqrt(d * J)
    grad_v <- matrix(0, nrow = J, ncol = d)
    grad_sigma2 <- 0

    for (j in seq_len(J)) {
      col_start <- (j - 1L) * d + 1L
      col_end <- j * d

      g_block <- grad_obj_tau[, col_start:col_end, drop = FALSE]
      vj <- matrix(v_mat[j, ], nrow = n_local, ncol = d, byrow = TRUE)
      delta <- x_train - vj
      sq_norm <- rowSums(delta * delta)
      k_xv <- exp(-sq_norm / (2 * sigma2_k))
      b <- grads - delta / sigma2_k

      a <- g_block * (scale_fac * k_xv)
      q <- rowSums(a * b)

      grad_v[j, ] <- colSums(delta * q) / sigma2_k + colSums(a) / sigma2_k

      c_coef <- sq_norm / (2 * sigma2_k^2)
      grad_sigma2 <- grad_sigma2 + sum(c_coef * q) + sum(a * delta) / (sigma2_k^2)
    }

    grad_log_sigma2 <- grad_sigma2 * sigma2_k
    grad_par <- c(as.numeric(t(grad_v)), grad_log_sigma2)

    list(value = fssd2 / denom, grad = grad_par)
  }

  objective_from_par <- function(par) {
    objective_and_grad_from_par(par)$value
  }

  grad_objective_from_par <- function(par) {
    objective_and_grad_from_par(par)$grad
  }

  opt <- stats::optim(
    par = par_init,
    fn = function(par) -objective_from_par(par),
    gr = function(par) -grad_objective_from_par(par),
    method = "L-BFGS-B",
    lower = c(rep(-Inf, J * d), log(1e-12)),
    upper = c(rep(Inf, J * d), log(1e6)),
    control = list(maxit = maxit)
  )

  par_opt <- opt$par
  v_opt <- matrix(par_opt[seq_len(J * d)], nrow = J, ncol = d, byrow = TRUE)
  sigma2_opt <- exp(par_opt[J * d + 1L])

  list(
    V_opt = v_opt,
    sigma2_opt = sigma2_opt,
    objective_opt = objective_from_par(par_opt),
    convergence = opt$convergence,
    message = opt$message
  )
}

#' Compute unbiased FSSD U-statistic on test split
#'
#' @param X_test Test samples (`n_te x d`).
#' @param grads_test Test scores (`n_te x d`).
#' @param V_opt Optimized test locations (`J x d`).
#' @param sigma2_opt Optimized kernel variance.
#'
#' @return List with `fssd_stat` and `tau_matrix`.
#' @export
compute_fssd_stat <- function(X_test, grads_test, V_opt, sigma2_opt) {
  x_test <- as_samples_matrix_fssd(X_test)
  grads <- as.matrix(grads_test)
  if (!is.numeric(grads) || nrow(grads) != nrow(x_test) || ncol(grads) != ncol(x_test)) {
    stop("grads_test must be a numeric matrix with the same shape as X_test")
  }

  tau <- compute_tau(x_test, grads, V_opt, sigma2_opt)
  fssd2 <- unbiased_fssd_stat(tau)

  list(fssd_stat = fssd2, tau_matrix = tau)
}

#' Simulate null distribution and compute FSSD p-value
#'
#' @param tau_matrix Matrix of test features (`n_te x dJ`).
#' @param fssd_stat Unbiased FSSD statistic from test split.
#' @param n_simulations Number of null simulations.
#'
#' @return List with p-value and null diagnostic values.
#' @export
fssd_pvalue <- function(tau_matrix, fssd_stat, n_simulations = 2000) {
  tau <- as.matrix(tau_matrix)
  if (!is.numeric(tau) || nrow(tau) < 2 || ncol(tau) < 1) {
    stop("tau_matrix must be a numeric matrix with at least two rows")
  }
  if (!is.numeric(fssd_stat) || length(fssd_stat) != 1 || !is.finite(fssd_stat)) {
    stop("fssd_stat must be a finite scalar")
  }
  if (!is.numeric(n_simulations) || length(n_simulations) != 1 ||
      !is.finite(n_simulations) || n_simulations <= 0) {
    stop("n_simulations must be a positive scalar")
  }
  n_simulations <- as.integer(n_simulations)

  n_te <- nrow(tau)
  sigma_q <- safe_second_moment(tau)

  eigvals <- eigen(sigma_q, symmetric = TRUE, only.values = TRUE)$values
  eigvals[eigvals < 0 & abs(eigvals) < 1e-10] <- 0
  eigvals[eigvals < 0] <- 0

  z <- matrix(stats::rnorm(n_simulations * length(eigvals)),
              nrow = n_simulations,
              ncol = length(eigvals))
  s_null <- as.numeric((z^2 - 1) %*% eigvals)

  test_score <- n_te * fssd_stat
  p_value <- mean(s_null > test_score)

  list(
    p_value = p_value,
    test_score = test_score,
    null_samples = s_null,
    eigenvalues = eigvals
  )
}

#' Compute tau feature matrix for FSSD
#'
#' @param X Samples (`n x d`).
#' @param grads Score matrix (`n x d`).
#' @param V Test locations (`J x d`).
#' @param sigma2_k Positive Gaussian kernel variance.
#'
#' @return Matrix of size `n x (dJ)`.
#' @export
compute_tau <- function(X, grads, V, sigma2_k) {
  x_mat <- as.matrix(X)
  grads_mat <- as.matrix(grads)
  v_mat <- as.matrix(V)

  if (!is.numeric(x_mat) || !is.numeric(grads_mat) || !is.numeric(v_mat)) {
    stop("X, grads, and V must be numeric matrices")
  }
  if (nrow(x_mat) != nrow(grads_mat) || ncol(x_mat) != ncol(grads_mat)) {
    stop("X and grads must have the same dimensions")
  }
  if (ncol(v_mat) != ncol(x_mat) || nrow(v_mat) < 1) {
    stop("V must be a numeric matrix with ncol(V) == ncol(X)")
  }

  if (!is.numeric(sigma2_k) || length(sigma2_k) != 1 || !is.finite(sigma2_k) || sigma2_k <= 0) {
    stop("sigma2_k must be a positive scalar")
  }

  n <- nrow(x_mat)
  d <- ncol(x_mat)
  J <- nrow(v_mat)
  scale_fac <- 1 / sqrt(d * J)

  tau <- matrix(0, nrow = n, ncol = d * J)

  for (j in seq_len(J)) {
    vj <- matrix(v_mat[j, ], nrow = n, ncol = d, byrow = TRUE)
    delta <- x_mat - vj
    sq_norm <- rowSums(delta * delta)
    k_xv <- exp(-sq_norm / (2 * sigma2_k))

    grad_k <- -(delta / sigma2_k) * k_xv
    xi <- grads_mat * k_xv + grad_k
    xi <- xi * scale_fac

    col_start <- (j - 1L) * d + 1L
    col_end <- j * d
    tau[, col_start:col_end] <- xi
  }

  tau
}

#' Coerce input samples to an `n x d` numeric matrix
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
  if (nrow(x) < 4) {
    stop("X must contain at least four samples")
  }

  x
}

#' Validate and coerce score output to an `n x d` matrix
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

#' Unbiased FSSD U-statistic via O(n) expansion
#' @keywords internal
unbiased_fssd_stat <- function(tau_matrix) {
  tau <- as.matrix(tau_matrix)
  n <- nrow(tau)
  if (n < 2) {
    stop("Need at least two rows to compute unbiased FSSD statistic")
  }

  sum_tau <- colSums(tau)
  norm_sum_sq <- sum(sum_tau * sum_tau)
  row_norm_sq <- rowSums(tau * tau)

  as.numeric((norm_sum_sq - sum(row_norm_sq)) / (n * (n - 1)))
}

#' Numerically stable covariance
#' @keywords internal
safe_covariance <- function(x) {
  x_mat <- as.matrix(x)
  p <- ncol(x_mat)
  if (p == 1) {
    v <- stats::var(x_mat[, 1])
    if (!is.finite(v) || v < 0) {
      v <- 0
    }
    return(matrix(v, nrow = 1, ncol = 1))
  }

  s <- stats::cov(x_mat)
  s[!is.finite(s)] <- 0
  s <- (s + t(s)) / 2
  s
}

#' Non-centered second moment matrix
#' @keywords internal
safe_second_moment <- function(x) {
  x_mat <- as.matrix(x)
  if (!is.numeric(x_mat) || nrow(x_mat) < 1 || ncol(x_mat) < 1) {
    stop("x must be a non-empty numeric matrix")
  }

  n <- nrow(x_mat)
  m2 <- crossprod(x_mat) / n
  m2[!is.finite(m2)] <- 0
  m2 <- (m2 + t(m2)) / 2
  m2
}

#' Median pairwise Euclidean distance heuristic
#' @keywords internal
median_pairwise_distance <- function(x, max_samples = 1000) {
  x_mat <- as.matrix(x)
  n <- nrow(x_mat)
  if (n < 2) {
    stop("Need at least two samples for median heuristic")
  }

  if (!is.numeric(max_samples) || length(max_samples) != 1 || !is.finite(max_samples) || max_samples < 2) {
    stop("max_samples must be a numeric scalar >= 2")
  }
  max_samples <- as.integer(max_samples)

  if (n > max_samples) {
    idx <- sample.int(n, max_samples)
    x_mat <- x_mat[idx, , drop = FALSE]
  }

  dmat <- as.matrix(stats::dist(x_mat, method = "euclidean"))
  vals <- dmat[upper.tri(dmat)]
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0) {
    return(1.0)
  }
  stats::median(vals)
}

#' Draw samples from a Gaussian with numerical stabilization
#' @keywords internal
sample_gaussian_rows <- function(n, mu, sigma) {
  d <- length(mu)
  sigma <- as.matrix(sigma)
  if (!all(dim(sigma) == c(d, d))) {
    stop("sigma must be d x d")
  }

  sigma <- (sigma + t(sigma)) / 2
  eig <- eigen(sigma, symmetric = TRUE)
  vals <- pmax(eig$values, 1e-10)
  root <- eig$vectors %*% diag(sqrt(vals), nrow = d)

  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  sweep(z %*% t(root), 2, mu, "+")
}
