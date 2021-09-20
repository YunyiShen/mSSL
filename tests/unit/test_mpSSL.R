library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)
library(mSSL)


source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/simu_data.R")
sourceCpp("./src/mpSSL.cpp")


set.seed(12345)
p <- 2
q <- 5
n <- 300
B <- as.matrix( rsparsematrix(p, q, 0.3, rand.x = function(n){runif(n,-2,2)}))

X <- matrix(rnorm(p * n), nrow = n, ncol = p)

graph <- g_model1(q)
Sigma <- graph$Sigma
Omega <- graph$Omega
mu <- runif(q,-2,2)

Y <- mprobit(X,B,mu,Sigma,Omega)

mpSSL_dpe_res <- mpSSL_dpe(X,Y,lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                    xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                    theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 5000,
                    eps = 1e-2,
                    s_max_condition = 10*nrow(X),
                    obj_counter_max = 5,
                    verbose = 1, n_rep = 200)

mpSSL_dcpe_res <- mpSSL_dcpe(X,Y,
                    lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                    xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                    theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 10000,
                    eps = 1e-6,
                    #s_max_condition = 10*nrow(X),
                    #obj_counter_max = 5,
                    verbose = 1, n_rep = 200)



CARlasso:::stein_loss(Omega, mpSSL_dcpe_res$Omega)
CARlasso:::stein_loss(Omega, mpSSL_dpe_res$Omega)



erB_dcpe <- error_B(GCSSL_dcpe_res$B, B)
erB_dpe <- error_B(mpSSL_dpe_res$B, B)

erB <- rbind(erB_dcpe,erB_dpe)
rownames(erB) <- c("cgSSL_dcpe","cgSSL_dpe")

erO_dcpe <- error_Omega(GCSSL_dcpe_res$Omega, Omega)
erO_dpe <- error_Omega(mpSSL_dpe_res$Omega, Omega)


erO <- rbind(erO_dcpe,erO_dpe)
rownames(erO) <- c("cgSSL_dcpe","cgSSL_dpe")

erB
erO


