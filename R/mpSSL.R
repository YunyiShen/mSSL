#' Multivariate probit Spike-and-Slab LASSO
#' @description Main posterior exploration algorithm for multivariate probit spike-and-slab LASSO.
#' @param Y response matrix
#' @param X design matrix
#' @param condexp bool, whether to do the fast conditional posterior exploration (dcpe), default is FALSE for doing the dynamic posterior exploration (dpe)
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
#' @return A list with dynamic exploration result, point estimates are in `$Omega` and `$B`.
#' 

mpSSL <- function(Y,X, condexp = FALSE, 
                  lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                  xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                  theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                  eta_hyper_params = c(1, ncol(Y)),
                  diag_penalty = FALSE,
                  max_iter = 500,
                  eps = 1e-3,
                  s_max_condition = 10*nrow(X),
                  obj_counter_max = 5,
                  verbose = FALSE){
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    diag_penalty <- 1 * diag_penalty

    if(length(unique(c(Y)))!=2){
        stop("Y must only have two unique values")
    }

    if(!all(sort(unique(c(Y)))!=c(0,1) )){
        warning("Y is not 0,1 coded, will convert larger value to 1")
        uniquevalues <- sort(unique(c(Y)))
        Y[Y==uniquevalues[1]] <- 0
        Y[Y==uniquevalues[2]] <- 1
    }


    if(nrow(X)!=nrow(Y)){
        stop("X and Y should have the number of rows.")
    }


    if(condexp){
        if(verbose){cat("start the dynamic *conditional* posterior exploration...\n")}
        verbose <- 1 * verbose
        res <- mpSSL_dcpe(X,Y, lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, verbose)
    }
    else {
        if(verbose){cat("start the dynamic posterior exploration...\n")}
        verbose <- 1 * verbose
        res <- mpSSL_dpe(X,Y, lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, s_max_condition, 
                    obj_counter_max,
                    verbose)
    }
    if(verbose==1){cat("done\n")}
    return(res)
}