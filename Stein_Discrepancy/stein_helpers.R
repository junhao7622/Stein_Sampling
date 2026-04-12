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

# Cleans up workspace by removing all variables.
cleanup <- function(){
  rm(list = ls())
}

# Score function for gamma distribution.
gamma_score <- function(x, shape, rate = 1, scale = 1 / rate){
  return((shape - 1) / x - 1 / scale)
}

# Function that can be used to retain legend of a plot.
get_legend <- function(myggplot){
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Function that rounds up given number to three significant digits.
custom_round <- function(x){
  round(x, 3)
}

# Median pairwise distance heuristic.
find_median_distance <- function(Z){
  if(is.data.frame(Z)){
    Z <- data.matrix(Z)
  }else{
    Z <- as.array(Z)
  }
  sample_size <- dim(Z)[1]
  num_dim <- dim(Z)[2]

  if(sample_size > 100){
    if(is.na(num_dim)){
      z_med <- Z[sample(sample_size, 100)]
    }else{
      z_med <- Z[sample(sample_size, 100), ]
    }
    sample_size <- 100
  }else{
    z_med <- Z
  }

  z_med_sq <- z_med * z_med
  if(is.na(dim(Z)[2])){
    gram_diag <- z_med_sq
  }else{
    gram_diag <- rowSums(z_med_sq)
  }

  q_mat <- rep_col(gram_diag, sample_size)
  r_mat <- rep_row(t(gram_diag), sample_size)

  dists <- q_mat + r_mat - 2 * z_med %*% t(z_med)
  dists[lower.tri(dists, diag = TRUE)] <- 0
  dists <- array(dists, dim = c(sample_size^2, 1))
  median_dist <- median(dists[dists > 0])

  return(median_dist)
}

# Ratio median heuristic.
ratio_median_heuristic <- function(Z, score_function){
  Z <- as.array(Z)
  sample_size <- dim(Z)[1]
  num_dim <- dim(Z)[2]

  if(sample_size > 100){
    if(is.na(num_dim)){
      z_med <- Z[sample(sample_size, 100)]
    }else{
      z_med <- Z[sample(sample_size, 100), ]
    }
    sample_size <- 100
    print("Sampled (Heuristic)")
  }else{
    z_med <- Z
    print("Original 100 dataset used (Heuristic)")
  }

  z_med_sq <- z_med * z_med
  if(is.na(dim(Z)[2])){
    gram_diag <- z_med_sq
  } else{
    gram_diag <- rowSums(z_med_sq)
  }

  q_mat <- rep_col(gram_diag, sample_size)
  r_mat <- rep_row(t(gram_diag), sample_size)

  dists <- q_mat + r_mat - 2 * z_med %*% t(z_med)
  dists[lower.tri(dists, diag = TRUE)] <- 0
  dists <- array(dists, dim = c(sample_size^2, 1))
  median_dist <- median(dists[dists > 0])

  if(is.na(num_dim)){
    z_med <- as.double(z_med)
  }
  score_x <- score_function(z_med)
  score_xy <- score_x %*% t(score_x)
  score_xy[lower.tri(score_xy, diag = TRUE)] <- 0
  score_xy <- array(score_xy, dim = c(sample_size^2, 1))
  median_score_xy <- median(score_xy[score_xy > 0])

  bandwidth <- (median_dist / median_score_xy)^(1 / 4)
  med_info <- list("h" = bandwidth, "median_dist" = median_dist, "median_sqxx" = median_score_xy)
  return(med_info)
}

# Legacy aliases (backward compatibility)
rep.row <- rep_row
rep.col <- rep_col
repmat <- rep_mat
getDim <- get_dim
custround <- custom_round
