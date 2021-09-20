#' Multivariate and Chain Graphical Spike-and-Slab LASSO for continuous data
#' @description Main posterior exploration algorithm for multivariate chain graphical spike-and-slab LASSO.
#' @param Y response matrix
#' @param X design matrix
#' @param cg bool, whether use the chain graphical parameterization, default true
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
#' @examples 
#' B <- diag(5)
#' graph <- g_modelAR1(5)
#' Omega <- graph$Omega
#' Sigma <- graph$Sigma
#' X <- matrix(rnorm(5*50),50,5)
#' Y <- rCARm(50,B,X,Sigma)
#' cgSSL_dpe_res <- mSSL(Y,X)
#' error_B(cgSSL_dpe_res$B,B)
#' error_Omega(cgSSL_dpe_res$Omega,Omega)

mSSL <- function(Y,X, cg = TRUE,
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
                  verbose = FALSE){
    Y <- as.matrix(Y)
    X <- as.matrix(X)
    diag_penalty <- 1 * diag_penalty
    if(nrow(X)!=nrow(Y)){
        stop("X and Y should have the number of rows.")
    }


    if(condexp){
        if(verbose){cat("start the dynamic *conditional* posterior exploration...\n")}
        verbose <- 1 * verbose
        if(cg){
            res <- cgSSL_dcpe(X,Y, lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, verbose)
        }
        else{
            res <- mSSL_dcpe(X,Y, lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, verbose)
        }
        
    }
    else {
        if(verbose){cat("start the dynamic posterior exploration...\n")}
        verbose <- 1 * verbose
        if(cg){
            res <- cgSSL_dpe(X,Y, lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, s_max_condition, 
                    obj_counter_max,
                    verbose)
        }
        else{
            res <- mSSL_dpe(X,Y, lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, s_max_condition, 
                    obj_counter_max,
                    verbose)
        }
        
    }
    if(verbose==1){cat("done\n")}
    return(res)
}


#' Graphical Spike-and-Slab LASSO
#' @description Main posterior exploration algorithm for graphical spike-and-slab LASSO, useful in general chain graph setting for the first layer
#' @param Y response matrix, no design in this case
#' @param xis hyperparameters to be explored by the algorithm, penalty on Omega
#' @param eta_hyper_params hyperparameter to be set, prior on spike weight of Omega
#' @param diag_penalty bool, whether to penalize the diagonal, default no
#' @param max_iter maximum iterations for the EM algorithm
#' @param eps tolerance for convergence
#' @param obj_counter_max , maximum number of couting the objective function
#' @param verbose bool, whether to print intermidate notes
#' @return A list with dynamic exploration result, point estimates are in `$Omega` and `$B`.
#' 

gSSL <- function(Y, 
                  xis = list(xi1 = 0.01 * nrow(Y), xi0 = seq(0.1 * nrow(Y), nrow(Y), length = 10)),
                  
                  eta_hyper_params = c(1, ncol(Y)),
                  diag_penalty = FALSE,
                  max_iter = 500,
                  eps = 1e-3,
                  
                  obj_counter_max = 5,
                  verbose = FALSE){
    Y <- as.matrix(Y)
    diag_penalty <- 1 * diag_penalty
    

    if(verbose){cat("start the dynamic posterior exploration...\n")}
    verbose <- 1 * verbose
    res <- gSSLcpp(Y, 
                    xis, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, 
                    obj_counter_max,
                    verbose)
    if(verbose==1){cat("done\n")}
    return(res)
}


#' Chain graphical Vector Autoregression Spike-and-Slab LASSO
#' @description Main posterior exploration algorithm for CVAR(1) spike-and-slab LASSO.
#' @param Y response matrix, a multivarate time series, row as observations
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

cgVARSSL <- function(Y, condexp = FALSE, 
                  lambdas = list(lambda1 = 1, lambda0 = seq(10, (nrow(Y)-1), length = 10)),
                    xis = list(xi1 = 0.01 * (nrow(Y)-1), xi0 = seq(0.1 * ((nrow(Y)-1)), (nrow(Y)-1), length = 10)),
                    theta_hyper_params = c(1, ncol(Y) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 10000,
                    eps = 1e-3,
                    s_max_condition = 10*(nrow(Y)-1),
                    obj_counter_max = 5,
                  verbose = FALSE){
    Y <- as.matrix(Y)
    diag_penalty <- 1 * diag_penalty


    if(condexp){
        if(verbose){cat("start the dynamic *conditional* posterior exploration...\n")}
        verbose <- 1 * verbose
        res <- cgVARSSL_dcpe(Y, lambdas, 
                    xis, theta_hyper_params, 
                    eta_hyper_params, 
                    diag_penalty, max_iter, 
                    eps, verbose)
    }
    else {
        if(verbose){cat("start the dynamic posterior exploration...\n")}
        verbose <- 1 * verbose
        res <- cgVARSSL_dpe(Y, lambdas, 
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

