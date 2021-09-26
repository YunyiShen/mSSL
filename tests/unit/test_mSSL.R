library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)
source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/plot_support.R")

sourceCpp("./src/mSSL.cpp",rebuild = TRUE)

set.seed(42)
p <- 2
q <- 5
n <- 100
B <- as.matrix( rsparsematrix(p, q, 0.2, rand.x = function(n){runif(n,-2,2)}))

X <- matrix(rnorm(p * n), nrow = n, ncol = p)
Y <- matrix(NA, n, q)

graph <- g_model1(q)

Sigma <- graph$Sigma
Omega <- graph$Omega

XB <- X %*% B

for(i in 1:n){
    Y[i,] <- mvrnorm(n = 1, XB[i,], Sigma)
}

mSSL_dpe_res <- mSSL_dpe(X,Y,lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                    xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                    theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 5000,
                    eps = 1e-2,
                    s_max_condition = 10*nrow(X),
                    obj_counter_max = 5,
                    verbose = 1)

mSSL_dpe_res <- mSSL_dcpe(X,Y,lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                    xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                    theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 5000,
                    eps = 1e-2,
                    verbose = 1)
