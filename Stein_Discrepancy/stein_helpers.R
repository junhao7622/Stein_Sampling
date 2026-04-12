# stein_helpers.R
# Shared helper functions for GMM and KSD implementations.

# Repeat a given vector row-wise.
rep_row <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

# Repeat a given vector column-wise.
rep_col <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

# Repeat a matrix in each dimension (Matlab-like repmat).
rep_mat <- function(x, a = 1, b = 1){
  rows <- dim(x)[1]
  cols <- dim(x)[2]
  if(is.null(cols)){
    cols <- 1
  }
  row_rep <- matrix(rep(t(x), a), ncol = cols, byrow = TRUE)
  new_x <- matrix(rep(row_rep, b), ncol = cols * b)
  return(new_x)
}

# Return sample size and dimension of x.
get_dim <- function(x){
  if(is.array(x)){
    n <- dim(x)[1]
    dimen <- dim(x)[2]
  }else{
    x <- array(x)
    n <- dim(x)[1]
    dimen <- 1
  }

  result <- list("n" = n, "dim" = dimen)
  return(result)
}

# Function that rounds up given number to three significant digits.
custom_round <- function(x){
  round(x, 3)
}

# Deterministic subsampling without leaking RNG state to callers.
sample_rows_seeded <- function(sample_size, take_n, seed = NULL) {
  if (is.null(seed)) {
    return(sample.int(sample_size, take_n))
  }

  if (requireNamespace("withr", quietly = TRUE)) {
    return(withr::with_seed(as.integer(seed), sample.int(sample_size, take_n)))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  sample.int(sample_size, take_n)
}

# Median pairwise squared-distance heuristic.
find_median_distance <- function(Z,
                                 max_samples = 2000,
                                 use_sampling = TRUE,
                                 seed = 123,
                                 distance_backend = c("auto", "base", "parallelDist")){
  if (!is.numeric(max_samples) || length(max_samples) != 1 || !is.finite(max_samples) || max_samples < 2) {
    stop("max_samples must be a numeric scalar >= 2")
  }
  max_samples <- as.integer(max_samples)
  if (!is.logical(use_sampling) || length(use_sampling) != 1 || is.na(use_sampling)) {
    stop("use_sampling must be TRUE or FALSE")
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 || !is.finite(seed))) {
    stop("seed must be NULL or a finite numeric scalar")
  }

  distance_backend <- match.arg(distance_backend)

  if (is.null(dim(Z))) {
    Z <- matrix(as.numeric(Z), ncol = 1)
  } else if (is.data.frame(Z)) {
    Z <- data.matrix(Z)
  } else {
    Z <- as.matrix(Z)
  }

  if (!is.numeric(Z) || nrow(Z) < 2) {
    stop("Z must be numeric with at least two rows")
  }

  sample_size <- nrow(Z)
  if (use_sampling && sample_size > max_samples) {
    idx <- sample_rows_seeded(sample_size, max_samples, seed = seed)
    z_med <- Z[idx, , drop = FALSE]
  } else {
    z_med <- Z
  }

  # Optional C++ backend through parallelDist for large problems.
  use_parallel_dist <- identical(distance_backend, "parallelDist") ||
    (identical(distance_backend, "auto") && requireNamespace("parallelDist", quietly = TRUE) && nrow(z_med) > 4000)

  if (use_parallel_dist && requireNamespace("parallelDist", quietly = TRUE)) {
    sq_dists <- as.numeric(parallelDist::parDist(z_med, method = "euclidean"))^2
  } else {
    sq_dists <- as.numeric(stats::dist(z_med, method = "euclidean"))^2
  }
  sq_dists <- sq_dists[is.finite(sq_dists)]
  if (length(sq_dists) == 0) {
    return(1)
  }

  med <- stats::median(sq_dists)
  if (!is.finite(med) || med < 0) {
    med <- 1
  }
  if (med == 0) {
    warning(
      "Median pairwise squared distance is zero (possible non-mixing/repeated states); using floor value 1e-5.",
      call. = FALSE
    )
    med <- 1e-5
  }

  med
}

# Legacy aliases (backward compatibility)
rep.row <- rep_row
rep.col <- rep_col
repmat <- rep_mat
getDim <- get_dim
custround <- custom_round
