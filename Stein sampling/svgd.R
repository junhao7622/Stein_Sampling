SVGD <- function() {
	self <- list()

	self$svgd_kernel <- function(theta, h = -1) {
		theta <- as.matrix(theta)

		sq_dist <- stats::dist(theta, method = "euclidean")
		pairwise_dists <- as.matrix(sq_dist)^2

		if (h < 0) {
			h <- stats::median(pairwise_dists)
			h <- sqrt(0.5 * h / log(nrow(theta) + 1))
		}

		# compute the rbf kernel
		Kxy <- exp(-pairwise_dists / h^2 / 2)

		dxkxy <- -Kxy %*% theta
		sumkxy <- rowSums(Kxy)
		for (i in seq_len(ncol(theta))) {
			dxkxy[, i] <- dxkxy[, i] + theta[, i] * sumkxy
		}
		dxkxy <- dxkxy / (h^2)

		list(Kxy = Kxy, dxkxy = dxkxy)
	}

	self$update <- function(x0,
													lnprob,
													n_iter = 1000,
													stepsize = 1e-3,
													bandwidth = -1,
													alpha = 0.9,
													debug = FALSE) {
		# Check input
		if (is.null(x0) || is.null(lnprob)) {
			stop("x0 or lnprob cannot be None!")
		}

		theta <- as.matrix(x0)

		# adagrad with momentum
		fudge_factor <- 1e-6
		historical_grad <- 0

		for (iter in seq_len(n_iter)) {
			if (debug && (iter %% 1000 == 0)) {
				message(paste0("iter ", iter))
			}

			lnpgrad <- lnprob(theta)

			# calculating the kernel matrix
			kernel_out <- self$svgd_kernel(theta, h = -1)
			kxy <- kernel_out$Kxy
			dxkxy <- kernel_out$dxkxy

			grad_theta <- (kxy %*% lnpgrad + dxkxy) / nrow(theta)

			# adagrad
			if (iter == 1) {
				historical_grad <- historical_grad + grad_theta^2
			} else {
				historical_grad <- alpha * historical_grad + (1 - alpha) * (grad_theta^2)
			}
			adj_grad <- grad_theta / (fudge_factor + sqrt(historical_grad))
			theta <- theta + stepsize * adj_grad
		}

		theta
	}

	self
}

