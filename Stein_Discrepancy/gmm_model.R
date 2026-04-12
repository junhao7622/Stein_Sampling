# gmm_model.R

if(file.exists("stein_helpers.R")){
  source("stein_helpers.R")
}

#' Returns a Gaussian Mixture Model
#'
#' @param nComp (scalar) : number of components
#' @param mu (d by k): mean of each component
#' @param sigma (d by d by k): covariance of each component
#' @param weights (1 by k) : mixing weight of each proportion (optional)
#' @param d : number of dimensions of vector (optional)
#'
#' @return model : A Gaussian Mixture Model generated from the given parameters
#'
#'
#' @export
#'
#' @examples
#' # Default 1-d gaussian mixture model
#' model <- gmm()
#'
#' # 1-d Gaussian mixture model with 3 components
#' model <- gmm(nComp = 3)
#'
#' # 3-d Gaussian mixture model with 3 components, with specified mu,sigma and weights
#' mu <- matrix(c(1,2,3,2,3,4,5,6,7),ncol=3)
#' sigma <- array(diag(3),c(3,3,3))
#' model <- gmm(nComp = 3, mu = mu, sigma=sigma, weights = c(0.2,0.4,0.4), d = 3)

gmm <- function(nComp=NULL, mu=NULL, sigma=NULL, weights=NULL, d=NULL){
  set.seed(0)
  # Generate a default Gaussian Mixture Model
  if(is.null(nComp)){
    d <- 1; k <- 5
    mu <- 10 * runif(k)
    mu <- mu - mean(mu)
    
    # Similate Matlab
    #mu <- c(0.9971,4.0301,-3.1720,-0.1498,-1.7055)
    
    sigma <- array(diag(d),c(d,d,k));
  }else{
    k <- nComp
  }
  
  # Case when mean is not specified
  if(is.null(mu)){
    if(is.null(d)) d <- 1
    
    mu <- 10 * replicate(k,runif(d))
    mu <- mu - mean(mu)
    
  }
  
  # Case when sigma is not specified
  if(is.null(sigma)){
    if(is.null(d)){
      d <- dim(mu)[1]
    }
    
    sigma <- array(diag(d),c(d,d,k));
  }
  
  #If sigma is one-dimensional
  if(is.null(dim(sigma))){
    sigma <- array(sigma, c(1,1,k))
  }
  
  # Handle the cases for the weights
  if(is.null(weights)){
    weights = rep(1,k) / k
  }else if(sum(weights < 0) > 1){
    stop('Non-positive weights')
  }else if(sum(weights) != 1){
    weights <- weights / sum(weights)
  }
  
  model <- list("nComp"=k, "mu" = mu, "sigma" = sigma,
                "weights" = weights, "d"=d)
  
  return(model)
}

#' Generates dataset from Gaussian Mixture Model
#' @note Requires library mvtnorm
#'
#' @param model : Gaussian Mixture Model defined by gmm()
#' @param n : number of samples desired
#'
#' @return data (n by d): Random dataset generated from given the Gaussian Mixture Model
#'
#'
#' @export
#'
#' @examples
#' #Generate 100 samples from default gaussian mixture model
#' model <- gmm()
#' X <- rgmm(model)
#'
#' #Generate 300 samples from 3-d gaussian mixture model
#' model <- gmm(d=3)
#' X <- rgmm(model,n=300)


rgmm <- function(model = NULL,n=100){
  if(is.null(model)){
    stop('Supply GMM Model')
  }else{
    k <- model$nComp
    mu <- model$mu
    sigma <- model$sigma
    weights <- model$weights
    d <- model$d
    
    components <- sample(1:k,prob=weights,size=n,replace=TRUE)
    
    # 1 Dimensional case
    if(d == 1){
      stdev <- rep(0,k)
      for(i in 1:k){
        stdev[i] <- sqrt(sigma[,,i])
      }
      data <- rnorm(n=n,mean=mu[components],sd=stdev[components])
    }else{
      # Multidimensional Case
      if (!requireNamespace("mvtnorm", quietly = TRUE)) {
        stop("mvtnorm needed for this demo to work. Please install it.",
             call. = FALSE)
      }
      data = NULL
      for(i in 1:k){
        count <- sum(components == i)
        single_cluster_data <- mvtnorm::rmvnorm(count, mean=mu[,i], sigma = sigma[,,i])
        data = rbind(data,single_cluster_data)
      }
    }
    
    return(data)
    # hcum <- h <- hist(data,breaks=30,plot=FALSE)
    # hcum$counts <- cumsum(hcum$counts)
  }
}


#' Returns a perturbed model of given GMM
#'
#' @param model : The base Gaussian Mixture Model
#'
#' @return perturbedModel : Perturbed model with added noise to the supplied GMM
#'
#'
#' @export
#'
#' @examples
#' #Add noise to default 1-d gaussian mixture model
#' model <- gmm()
#' noisymodel <- perturbgmm(model)

perturbgmm <- function(model = NULL){
  if(is.null(model)){
    stop('Supply GMM Model')
  }
  set.seed(0)
  
  k <- model$nComp
  d <- model$d
  noise <- replicate(k,rnorm(d))
  perturbed_mu <- model$mu + noise
  perturbed_model <- gmm(k,perturbed_mu,model$sigma, model$weights,d)
  
  return(perturbed_model)
}

#' Calculates the posterior probability for a given dataset for a GMM
#'
#' @param model : The Gaussian Mixture Model
#' @param X (n by d): The dataset of interest,
#' where n is the number of samples and d is the dimension
#'
#' @return P (n by k) : The posterior probabilty of each dataset belonging to each of the k component
#'
#'
#' @export
#'
#' @examples
#' # compute posterior probability for a default 1-d gaussian mixture model
#' # and dataset generated from it
#' model <- gmm()
#' X <- rgmm(model)
#' p <- posteriorgmm(model=model, X=X)

posteriorgmm <- function(model=NULL, X=NULL){
  if(is.null(model) || is.null(X)){
    stop('Supply Model and Data')
  }else{
    
    n <- dim(X)[1];
    if(is.null(n)){
      n <- length(X)
    }
    d <- model$d;
    k <- model$nComp
    mu <- model$mu
    sigma <- model$sigma
    weights <- model$weights
    
    posterior_prob <- NULL
    # If one-dimensional
    if(d==1){
      stdev <- rep(0,k)
      for(i in 1:k){
        stdev <- sqrt(sigma[,,i])
        single_model_prob <- dnorm(X, mean=mu[i],sd <- stdev)
        posterior_prob <- cbind(posterior_prob, single_model_prob)
      }
    }else{
      if (!requireNamespace("mvtnorm", quietly = TRUE)) {
        stop("mvtnorm needed for this demo to work. Please install it.",
             call. = FALSE)
      }
      for(i in 1:k){
        single_model_prob <- mvtnorm::dmvnorm(X,mean=mu[,i],sigma=sigma[,,i])
        posterior_prob <- cbind(posterior_prob, single_model_prob)
      }
    }
    
    posterior_prob <- posterior_prob * rep_row(weights,n)
    sum_prob <- rowSums(posterior_prob)
    posterior_prob <- posterior_prob / replicate(k,sum_prob)
    posterior_prob[is.nan(posterior_prob)] <- 0
  }
  
  return(posterior_prob)
}

#' Calculates the likelihood for a given dataset for a GMM
#'
#' @param model : The Gaussian Mixture Model
#' @param X (n by d): The dataset of interest,
#' where n is the number of samples and d is the dimension
#'
#' @return P (n by k) : The likelihood of each dataset belonging to each of the k component
#'
#'
#' @export
#'
#' @examples
#' # compute likelihood for a default 1-d gaussian mixture model
#' # and dataset generated from it
#' model <- gmm()
#' X <- rgmm(model)
#' p <- likelihoodgmm(model=model, X=X)
likelihoodgmm <- function(model=NULL, X=NULL){
  if(is.null(model) || is.null(X)){
    stop('Supply Model and Data')
  }else{
    n <- dim(X)[1];
    if(is.null(n)){
      n <- length(X)
    }
    d <- model$d;
    k <- model$nComp
    mu <- model$mu
    sigma <- model$sigma
    weights <- model$weights
    
    likelihood_prob <- NULL
    # If one-dimensional
    if(d==1){
      stdev <- rep(0,k)
      for(i in 1:k){
        stdev <- sqrt(sigma[,,i])
        single_model_prob <- dnorm(X, mean=mu[i],sd <- stdev)
        likelihood_prob <- cbind(likelihood_prob, single_model_prob)
      }
    }else{
      if (!requireNamespace("mvtnorm", quietly = TRUE)) {
        stop("mvtnorm needed for this demo to work. Please install it.",
             call. = FALSE)
      }
      
      for(i in 1:k){
        single_model_prob <- mvtnorm::dmvnorm(X,mean=mu[,i],sigma=sigma[,,i])
        likelihood_prob <- cbind(likelihood_prob, single_model_prob)
      }
    }
    
    likelihood_prob <- likelihood_prob * rep_row(weights,n)
    sum_prob <- rowSums(likelihood_prob)
  }
  
  return(sum_prob)
}


#' Score function for given GMM : calculates score function dlogp(x)/dx
#' for a given Gaussian Mixture Model
#'
#' @param model : The Gaussian Mixture Model
#' @param X (n by d): The dataset of interest,
#' where n is the number of samples and d is the dimension
#'
#' @return y : The score computed by the given function
#'
#'
#' @export
#'
#' @examples
#' # Compute score for a given gaussianmixture model and dataset
#' model <- gmm()
#' X <- rgmm(model)
#' score <- scorefunctiongmm(model=model, X=X)
scorefunctiongmm <- function(model=NULL, X=NULL){
  if(is.null(model) || is.null(X)){
    stop('Supply Model and Data')
  }else{
    P <- posteriorgmm(model, X)
    d <- model$d
    if(d == 1){
      n <- length(X)
    }else{
      n <- dim(X)[1]
    }
    
    score <- 0
    
    for(component_idx in 1:model$nComp){
      sigma <- model$sigma[,,component_idx]
      if(d == 1){
        score <- score + replicate(d,P[,component_idx]) * t(qr.solve(sigma,t(-X + rep_row(model$mu[component_idx],n))))
      }else{
        score <- score + replicate(d,P[,component_idx]) * t(qr.solve(sigma,t(-X + rep_row(model$mu[,component_idx],n))))
      }
      
    }
  }
  
  return(score)
}


#' Plots histogram for 1-d GMM given the dataset
#'
#' @param data (n by 1): The dataset of interest,
#' where n is the number of samples.
#' @param mu : True mean of the GMM (optional)
#'
#'
#' @export
#'
#' @examples
#' # Plot pdf histogram for a given dataset
#' model <- gmm()
#' X <- rgmm(model)
#' plotgmm(data=X)
#'
#' # Plot pdf histogram for a given dataset, with lines that indicate the mean
#' model <- gmm()
#' mu <- model$mu
#' X <- rgmm(model)
#' plotgmm(data=X, mu=mu)
plotgmm <- function(data, mu = NULL){
  hcum <- h <- hist(data,breaks=30,plot=FALSE)
  hcum$counts <- cumsum(hcum$counts)
  
  ##----- Plot only pdf----------------
  par(mar=c(1,1,1,1))
  plot(h, xlim = c(-20,20),col="grey")
  d <- density(data)
  lines(x = d$x, y = d$y  * length(data) * diff(h$breaks)[1], lwd = 2)
  if(!is.null(mu)){
    abline(v = mu,col = "royalblue",lwd = 2)
  }
}


#' Returns a score evaluator closure for a fixed GMM
#'
#' This helper adapts a generated GMM model to interfaces that require
#' `grad_log_prob` with signature `f(X)`.
#'
#' @param model : Gaussian Mixture Model defined by gmm()
#'
#' @return A function that takes only `X` and returns gradients of log-density.
#' If `X` is a vector (no `dim`) and `model$d == 1`, it is treated as a batch
#' of 1-D samples and the returned gradient is a vector of the same length.
#' If `X` is a vector and `model$d > 1`, it is treated as a single sample in
#' `d` dimensions and the returned gradient is a vector of length `d`.
#' If `X` is a matrix, the returned gradient is a matrix.
#'
#' @export
#'
#' @examples
#' model <- gmm()
#' grad_log_prob <- get_score_evaluator(model)
#' X <- rgmm(model)
#' G <- grad_log_prob(X)
get_score_evaluator <- function(model){
  function(X){
    if(is.null(dim(X))){
      if(model$d == 1){
        X_mat <- matrix(X, ncol = 1)
      }else{
        X_mat <- matrix(X, nrow = 1)
      }
      grad <- scorefunctiongmm(model = model, X = X_mat)
      return(as.vector(grad))
    }

    scorefunctiongmm(model = model, X = X)
  }
}

# Snake_case aliases
gaussian_mixture_model <- gmm
sample_gmm <- rgmm
perturb_gmm <- perturbgmm
posterior_gmm <- posteriorgmm
likelihood_gmm <- likelihoodgmm
score_function_gmm <- scorefunctiongmm
plot_gmm <- plotgmm