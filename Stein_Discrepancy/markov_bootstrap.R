# markov_bootstrap.R

library(stats)  # Provides acf().

simulate_ar <- function(nPeriod, nPath, beta) {
  noise <- matrix(rnorm(nPeriod * nPath), nrow = nPeriod, ncol = nPath)
  sims  <- matrix(0,                      nrow = nPeriod, ncol = nPath)
  sims[1, ] <- noise[1, ]
  sqrt_beta  <- sqrt(1 - beta^2)
  for (period in 2:nPeriod) {
    sims[period, ] <- beta * sims[period - 1, ] + sqrt_beta * noise[period, ]
  }
  return(sims)
}

simulatepm <- function(N, p_change) {
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

# plot main
if (sys.nframe() == 0) {
  w <- simulatepm(1000, 0.02)
  print(acf(w, lag.max = 3, plot = FALSE)$acf)
  print(as.numeric(t(w) %*% w) / 10000)

  plot(w, type = "l", main = "simulatepm")
}
