#' Multivariate probit model
#' simulating from a multivariate probit model
#' @param X design matrix
#' @param B regression coefficients
#' @param mu mean vector
#' @param Sigma covariance
#' @param Omega precision matrix
#' @return a matrix of binary variables sampled from the multivariate probit with mean-covariance parameterization 
mprobit <- function(X,B,mu,Sigma, Omega=solve(Sigma)){
    unitdiag(Sigma,Omega)
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    1.*(XB+E>=0)

}

#' Multivariate probit model with chain graph parameterization 
#' simulating from a multivariate probit model
#' @param X design matrix
#' @param B regression coefficients
#' @param mu mean vector
#' @param Sigma covariance
#' @param Omega precision matrix
#' @return a matrix of binary variables sampled from the multivariate probit with chain graph parameterization 
cgprobit <- function(X,B,mu,Sigma, Omega=solve(Sigma)){
  unitdiag(Omega,Sigma)
  XB <- X %*% B
  E <- mvrnorm(nrow(X),mu,Omega)
  res <- XB+E
  res <- res %*% Sigma
  1.*(res>=0)
}


#' Sampling from CAR, simplest Gaussian chain graph model
#' @description Take random samples following CAR model, i.e. N(XBSigma,Sigma)
#' @param n number of samples
#' @param B regression coefficients
#' @param X Design matrix, should have more than n rows
#' @param Sigma Covariance matrix
#' @return a matrix of n rows

rCARm <- function(n,B,X,Sigma){
    if(nrow(X)<n){
        stop("Please provide enough rows in design matrix")
    }

    if(nrow(B)!=ncol(X)){
        stop("Dimension mismatch of B and X")
    }
    if(ncol(B)!=ncol(Sigma)){
        stop("Dimension mismatch of B and Sigma")
    }
    res <- matrix(NA, nrow = nrow(X), ncol = ncol(Sigma))
    XB <- X %*% B
    for(i in 1:n){
        res[i,] <- MASS::mvrnorm(1,Sigma %*% (XB[i,]), Sigma)
    }
    return(res)
}



#' Sampling from Conditional Vector Autoregression (1) model, another simple Gaussian chain graph model
#' @description Take random samples following VCAR(1) model, i.e. Y_t+1~N(Y_tBSigma,Sigma)
#' @param n number of samples
#' @param B regression coefficients
#' @param Sigma Covariance matrix
#' @return a matrix of n rows

rCVARm <- function(n, B,Sigma){
    if(nrow(B)!=ncol(B) | nrow(B)!=ncol(Sigma)){
        stop("Dimension mismatch of B and Sigma")
    }
    q <- nrow(Sigma)
    res <- matrix(NA,n,q)
    res[1,] <- MASS::mvrnorm(1,rep(0,q),Sigma)
    for(i in 2:n){
        XB <- res[i-1,] %*% B
        res[i,] <- MASS::mvrnorm(1,XB %*% Sigma,Sigma)
    }
    return(res)
}

#' Covariace and precision matrix of AR(1) model
#' generating precision matrix and covariance of an AR(1) model
#' @param q number of nodes
#' @param rho the covariance of neighborhood
#' @param tol tolerence of 0 in precision
#' @return a list with Omega (precision) and Sigma (covariance)
g_modelAR1 <- function(q, rho=.7, tol = 1e-10){
    temp <- matrix(rep(1:q,q),ncol = q)
    Sigma <- rho ^ (abs(temp-t(temp)))
    Omega <- solve(Sigma)
    Omega <- Omega * (abs(Omega)>tol)
    return(list(Sigma = Sigma, Omega = Omega))
}

rSTARlogfloor <- function(X,B,mu,Sigma){
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    res <- XB+E
    res <- floor(exp(res))
    return(res)
}

rSTARidfloor0 <- function(X,B,mu,Sigma){
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    res <- XB+E
    res <- floor(res)
    return(res * (res>=0))
}


rSTARidfloor <- function(X,B,mu,Sigma){
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    res <- XB+E
    res <- floor(res)
    return(res)
}

#' simulating a STAR model with latent normal being mean-covariance parameterization
#' @description Take random samples following STAR and mean-covariance parameterization of latent normal, using link function given 
#' @param X design matrix
#' @param B regression coefficients
#' @param mu mean vector
#' @param Sigma covariance
#' @param link a function convert latent normal to integers, default is the log-floor link
#' @return a matrix of STAR sample
rSTAR <- function(X,B,mu, Sigma,link = function(x){floor(exp(x))}){
    if(nrow(B)!=ncol(X)){
        stop("Dimension mismatch of B and X")
    }
    if(ncol(B)!=ncol(Sigma)){
        stop("Dimension mismatch of B and Sigma")
    }
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    res <- XB+E
    return(link(res))
}


#' simulating a STAR model with latent normal being chain graph parameterization
#' @description Take random samples following STAR and chain graph parameterization of latent normal, using link function given 
#' @param X design matrix
#' @param B conditional regression coefficients
#' @param mu conditional mean vector
#' @param Sigma covariance
#' @param link a function convert latent normal to integers, default is the log-floor link
#' @return a matrix of STAR sample
rcgSTAR <- function(X,B,mu, Sigma,link = function(x){floor(exp(x))}){
    if(nrow(B)!=ncol(X)){
        stop("Dimension mismatch of B and X")
    }
    if(ncol(B)!=ncol(Sigma)){
        stop("Dimension mismatch of B and Sigma")
    }
    Omega <- solve(Sigma)
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Omega)
    res <- XB+E
    res <- res %*% Sigma
    return(link(res))
}
