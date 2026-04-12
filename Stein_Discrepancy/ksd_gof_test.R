if (file.exists("ksd_u_test.R")) {
  source("ksd_u_test.R")
} else {
  stop("Required file 'ksd_u_test.R' not found in working directory.")
}

if (file.exists("ksd_v_test.R")) {
  source("ksd_v_test.R")
} else {
  stop("Required file 'ksd_v_test.R' not found in working directory.")
}

if (file.exists("fssd_test.R")) {
  source("fssd_test.R")
} else {
  stop("Required file 'fssd_test.R' not found in working directory.")
}

#' Bootstrap Stein Goodness-of-Fit Test (U, V, or FSSD)
#'
#' High-level wrapper for Stein discrepancy goodness-of-fit tests.
#'
#' For `stat_type = "U"`, this follows Algorithm 1 from Liu et al. (2016):
#' weighted multinomial bootstrap and uses the returned `htest` p-value
#' from `ksd_u_test`.
#'
#' For `stat_type = "V"`, this uses `ksd_v_test` with
#' `test_type = "aggregated"` and reads the p-value directly from the
#' returned `htest` object.
#'
#' For `stat_type = "FSSD"`, this uses `fssd_test`, a linear-time Finite
#' Set Stein Discrepancy test with a chi-squared reference p-value.
#'
#' Hypotheses:
#' H0: the sample set {x_i} is drawn from q (target distribution).
#' H1: the sample set {x_i} is not drawn from q.
#'
#' @param X Numeric vector or matrix of samples. A vector is treated as `n x 1`.
#' @param score_function Function that evaluates target distribution
#' @param m Number of bootstrap replications. Defaults to `1000`.
#' @param alpha Significance level in (0, 1). Defaults to `0.05`.
#' @param stat_type Statistic type. `"U"`, `"V"`, or `"FSSD"`. Defaults to `"U"`.
#' @param ... Additional arguments passed to the selected backend test.
#'
#' @return A list with entries:
#' `decision`: "Reject H0" or "Fail to reject H0".
#' `p_value`: p-value from the selected test backend.
#' `statistic`: Observed KSD statistic for the chosen `stat_type`.
#' `alpha`: Significance level used.
#' `stat_type`: Statistic type used (`"U"`, `"V"`, or `"FSSD"`).
#'
#' @examples
#' # model <- gmm(nComp = 3, d = 2)
#' # X <- rgmm(model, n = 80)
#' # grad_log_prob <- get_score_evaluator(model)
#' # ksd_bootstrap_gof_test(X, grad_log_prob, m = 500, alpha = 0.05, stat_type = "U")
#' # ksd_bootstrap_gof_test(X, grad_log_prob, m = 500, alpha = 0.05, stat_type = "V")
#' # ksd_bootstrap_gof_test(X, grad_log_prob, stat_type = "FSSD",
#' #                        num_random_frequencies = 5, scaling = c(1, 10))
#'
#' @export
ksd_bootstrap_gof_test <- function(X, score_function, m = 1000, alpha = 0.05,
                                   stat_type = c("U", "V", "FSSD"), ...) {
  stat_type <- match.arg(stat_type)

  if (!is.function(score_function)) {
    stop("score_function must be a function")
  }
  if (!is.numeric(m) || length(m) != 1 || !is.finite(m) || m <= 0) {
    stop("m must be a positive scalar")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0, 1)")
  }

  m <- as.integer(m)

  if (stat_type == "U") {
    test_result <- ksd_u_test(
      X = X,
      score_function = score_function,
      boot_method = "weighted",
      nboot = m,
      ...
    )
  } else if (stat_type == "V") {
    test_result <- ksd_v_test(
      X = X,
      score_function = score_function,
      test_type = "aggregated",
      nboot = m,
      ...
    )
  } else {
    test_result <- fssd_test(
      X = X,
      score_function = score_function,
      ...
    )
  }

  observed_stat <- as.numeric(test_result$statistic[[1]])
  p_value <- as.numeric(test_result$p.value)

  decision <- if (p_value < alpha) {
    "Reject H0"
  } else {
    "Fail to reject H0"
  }

  list(
    decision = decision,
    p_value = p_value,
    statistic = observed_stat,
    alpha = alpha,
    stat_type = stat_type
  )
}
