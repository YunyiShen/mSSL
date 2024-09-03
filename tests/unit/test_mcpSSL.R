library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)



source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/simu_data.R")
sourceCpp("./src/mSSL.cpp",rebuild = T)



set.seed(12)
p <- 2
q <- 5
n <- 200
B <- as.matrix( rsparsematrix(p, q, 0.3, rand.x = function(n){runif(n,-1,1)}))

X <- matrix(rnorm(p * n), nrow = n, ncol = p)

graph <- g_model1(q)
Sigma <- graph$Sigma
Omega <- graph$Omega
mu <- runif(q,-1,1)

# first two are binary
#normalizing <- sqrt(Omega[1,1])
#Omega[1,] <- Omega[1,]/normalizing
#Omega[,1] <- Omega[,1]/normalizing


#normalizing <- sqrt(Omega[2,2])
#Omega[2,] <- Omega[2,]/normalizing
#Omega[,2] <- Omega[,2]/normalizing

#Sigma <- solve(Omega)
set.seed(42)
Y <- Y_latent <- X%*%B + rep(1, n) %*% t(mu) + mvrnorm(n, rep(0, q), Sigma = Sigma)
Y[,1:2] <- 1. * (Y[,1:2] >= 0)


#Y <- mprobit(X,B,mu,Sigma,Omega)

mpcSSL_dpe_res <- mpcSSL_dpe(X,Y,binidxend = 1,lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                    xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                    theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 500,
                    eps = 1e-3,
                    s_max_condition = 10 * nrow(X),
                    obj_counter_max = 5,
                    verbose = 1, nrep = 1000, nskp = 1)


erB_dpe <- error_B(mpcSSL_dpe_res$B, B)

erO_dpe <- error_Omega(mpcSSL_dpe_res$Omega, Omega)



