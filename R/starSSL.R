
#' the logfloor transformation to integer values
#' calculate the latent normal's bound under log-floor ransformation, i.e. y=floor(exp(z*)), z*~N
#' @param Y the observed data
#' @param upper whether to calculate upper bound
#' @return a matrix of lower or upper bound of the latent normal
logfloorlink <- function(Y, upper = TRUE){
    if(upper){
        return(log(Y+1))
    }
    return(log(Y))
}

#' identity floor 0 truncation link
#' calculate the latent normal's bound under identity floor 0 truncation link ransformation, i.e. y=floor(z*), y=y*(y>=0), z*~N
#' @param Y the observed data
#' @param upper whether to calculate upper bound
#' @return a matrix of lower or upper bound of the latent normal
idfloor0link <- function(Y,upper = TRUE){
    if(min(Y)<0) stop("Y should be non-negative, consider using idfloor link")
    if(upper){
        return(Y+1)
    }
    res <- Y
    res[Y==0] <- -Inf
    return(res)
}

#' identity floor link
#' calculate the latent normal's bound under identity floor link ransformation, i.e. y=floor(z*), y=y*(y>=0), z*~N
#' @param Y the observed data
#' @param upper whether to calculate upper bound
#' @return a matrix of lower or upper bound of the latent normal
idfloorlink <- function(Y,upper = TRUE){
    if(upper){
        return(Y+1)
    }
    return(Y)
}


#' Multivariate STAR Spike-and-Slab LASSO with fixed link function
#' @description Main posterior exploration algorithm for multivariate STAR spike-and-slab LASSO with fixed link function and floor rouding.
#' @param Y response matrix
#' @param X design matrix
#' @param link the function calculate the bounds of the latent normal h^-1(g(y*)), default is floor rounding with log link, should be a function with two arguments, the data and a boolean argument, if true we calculate the upper bound and lower bound otherwise.
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
                  nskp = 1){
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    diag_penalty <- 1 * diag_penalty


    if(nrow(X)!=nrow(Y)){
        stop("X and Y should have the number of rows.")
    }

    upper <- link(Y,TRUE)
    lower <- link(Y,FALSE)


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
    return(res)

}