
#' the logfloor transformation to integer values
#' calculate the latent normal's bound under log-floor ransformation, i.e. y=floor(exp(z*)), z*~N
#' @param Y the observed data
#' @return a list with both lower or upper bound of the latent normal
logfloorlink <- function(Y){
    return(list(lower = log(Y),upper = log(Y+1), invtransforming = function(z){floor(exp(z))}))
}

#' identity floor 0 truncation link
#' calculate the latent normal's bound under identity floor 0 truncation link ransformation, i.e. y=floor(z*), y=y*(y>=0), z*~N
#' @param Y the observed data
#' @return a list of lower and upper bound of the latent normal
idfloor0link <- function(Y){
    if(min(Y)<0) stop("Y should be non-negative, consider using idfloor link")
    upper <- (Y+1)
    lower <- Y
    lower[Y==0] <- -Inf
    return(
        list(lower = lower, 
            upper = upper, 
            invtransforming = function(t){max(floor(t),0)}
            )
        )
}

#' square root floor link truncated at 0
#' calculate the lower and upper bound when using a square root floor link
#' @param Y the data
#' @return a list with lower and upper bound and inverse transformation function 
sqrtfloor0link <- function(Y){
    lower <- sqrt(Y)
    lower[Y==0] <- -Inf
    upper <- sqrt(Y+1)

    return(list(lower = lower, 
                upper = upper, 
                invtransforming = function(w){floor(w^2)}))
}

#' identity floor link
#' calculate the latent normal's bound under identity floor link ransformation, i.e. y=floor(z*), y=y*(y>=0), z*~N
#' @param Y the observed data
#' @return a list of lower and upper bound of the latent normal
idfloorlink <- function(Y){
    return(list(lower = Y, 
                upper = Y+1, 
                invtransforming = function(w){w}))
}




#' Empirical link b Kowal & Wu 2021 and a floor rounding
#' calculate the lower and upper bound by under Kowal & Wu 2021's transformation deterimination procedure
#' @param Y the raw counting data
#' @param distribution which distribution to use?  must be one of
#' \itemize{
#' \item "np" (empirical CDF)
#' \item "pois" (moment-matched marginal Poisson CDF)
#' \item "neg-bin" (moment-matched marginal Negative Binomial CDF)
#' }
#' @param grid the grid to calculate the inverse transformation
#' @return A list of lower and upper bounds, as well as the empirical paramters used in that transformation, to be used for prediction.
KowalWufloorlink <- function(Y, distribution = "np", grid = NULL){
    p <- ncol(Y)
    trans_funs <- lapply(1:p, function(i, Y, distribution){
        g_cdf(Y[,i], distribution)
    }, Y, distribution)

    lower <- lapply(1:p, function(i,Y,trans_funs){
        trans_funs[[i]](Y[,i])
    },Y,trans_funs)
    lower <- Reduce(cbind, lower)
    upper <- lapply(1:p, function(i,Y,trans_funs){
        trans_funs[[i]](Y[,i] + 1)
    },Y,trans_funs)
    upper <- Reduce(cbind, upper)
    if(is.null(grid)){
        grid <- seq(min(5*min(lower),-1),max(5*max(upper),1),by = .1)
    }
    invtransforming <- lapply(trans_funs, g_inv_approx, grid)

    return(list(lower = lower, upper = upper, invtransforming = invtransforming))

}

## empirical transformation by Kowal & Wu, from rSTAR
g_cdf <- function(y, distribution = "np") {

  # Check: does the distribution make sense?
  distribution = tolower(distribution);
  if(!is.element(distribution, c("np", "pois", "neg-bin", "box-cox")))
    stop("The distribution must be one of 'np', 'pois', or 'neg-bin'")

  # Number of observations:
  n <- length(y)

  # Moments of the raw counts:
  mu_y <- mean(y); sigma_y = sd(y)

  # CDFs:
  if(distribution == 'np') {
    # (Scaled) empirical CDF:
    F_y <- function(t) n/(n+1)*ecdf(y)(t)
  }
  if(distribution == 'pois'){
    # Poisson CDF with moment-matched parameters:
    F_y = function(t) ppois(t,
                            lambda = mu_y)
  }
  if(distribution == 'neg-bin') {
    # Negative-binomial CDF with moment-matched parameters:
    if(mu_y >= sigma_y^2){
      # Check: underdispersion is incompatible with Negative-Binomial
      warning("'neg-bin' not recommended for underdispersed data")

      # Force sigma_y^2 > mu_y:
      sigma_y = 1.1*sqrt(abs(mu_y))
    }
    F_y = function(t) pnbinom(t,
                              size = mu_y^2/(sigma_y^2 - mu_y),
                              prob = mu_y/sigma_y^2)
  }

  # Input points for smoothing:
  t0 <- sort(unique(y[y!=0]))

  # Initial transformation:
  g0 <- mu_y + sigma_y*qnorm(F_y(t0-1))

  # Make sure we have only finite values of g0 (infinite values occur for F_y = 0 or F_y = 1)
  t0 <- t0[which(is.finite(g0))]; g0 = g0[which(is.finite(g0))]

  # Return the smoothed (monotone) transformation:
  splinefun(t0, g0, method = 'monoH.FC')
}

## Approximate inverse transformation by Kowal & Wu, from rSTAR
g_inv_approx <- function(g, t_grid) {

  # Evaluate g() on the grid:
  g_grid <- g(t_grid)

  # Approximate inverse function:
  function(s) {
    sapply(s, function(si)
      t_grid[which.min(abs(si - g_grid))])
  }
}

#' Multivariate STAR Spike-and-Slab LASSO with fixed link function
#' @description Main posterior exploration algorithm for multivariate STAR spike-and-slab LASSO with fixed link function and floor rouding.
#' @param Y response matrix
#' @param X design matrix
#' @param link the function calculate the bounds of the latent normal h^-1(g(y*)), default is floor rounding with log link, should be a function takes data and return a list with at least components lower and upper, with additional components for prediction.
#' @param cg bool whether to use chain graph parameterization, default FALSE
#' @param condexp bool, whether to do the fast conditional posterior exploration (dcpe), default is FALSE for doing the dynamic posterior exploration (dpe), not very stable when cg=TRUE
#' @param lambdas hyperparameters to be explored by the algorithm, penalty on B
#' @param xis hyperparameters to be explored by the algorithm, penalty on Omega
#' @param theta_hyper_params hyperparameter to be set, prior on spike weight of B
#' @param eta_hyper_params hyperparameter to be set, prior on spike weight of Omega
#' @param diag_penalty bool, whether to penalize the diagonal, default no
#' @param max_iter maximum iterations for the EM algorithm
#' @param eps tolerance for convergence
#' @param s_max_condition only used in dpe, maximum tolerance for the condition number
#' @param obj_counter_max only used in dpe, maximum number of couting the objective function
#' @param verbose bool, whether to print intermidate notes
#' @param nrep number of sample to take during E step of the latent normal
#' @param nskp one point taken per nskp samples during E step of the latent normal
#' @param ... other parameters for link function

#' @return A list with dynamic exploration result, point estimates are in `$Omega` and `$B`.
#' 
fstarmSSL <- function(Y,X,link = logfloorlink,
                  cg = FALSE, 
                  condexp = FALSE, 
                  lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                  xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                  theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                  eta_hyper_params = c(1, ncol(Y)),
                  diag_penalty = FALSE,
                  max_iter = 500,
                  eps = 1e-3,
                  s_max_condition = 10*nrow(X),
                  obj_counter_max = 5,
                  verbose = FALSE,
                  nrep = 200,
                  nskp = 1, ...){
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    diag_penalty <- 1 * diag_penalty


    if(nrow(X)!=nrow(Y)){
        stop("X and Y should have the number of rows.")
    }

    the_link <- link(Y,...)
    upper <- the_link$upper
    lower <- the_link$lower


    # mean-cov parameterization 
    if(condexp){
        if(verbose){cat("start the dynamic *conditional* posterior exploration...\n")}
        verbose <- 1 * verbose
        if(cg){
            res <- cgfstarSSL_dcpe(X,lower, upper,lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, verbose,nrep,nskp)
        }
        else{
            res <- fstarSSL_dcpe(X,lower, upper,lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, verbose,nrep,nskp)
        }
    }
    else {
        if(verbose){cat("start the dynamic posterior exploration...\n")}
        verbose <- 1 * verbose
        if(cg){
            res <- cgfstarSSL_dpe(X,lower,upper,lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, s_max_condition, 
                    obj_counter_max,
                    verbose,nrep,nskp)
        }
        else{
            res <- fstarSSL_dpe(X,lower,upper,lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, s_max_condition, 
                    obj_counter_max,
                    verbose,nrep,nskp)
        }
        
    }
    if(verbose==1){cat("done\n")}
    res$link_res <- the_link
    return(res)

}