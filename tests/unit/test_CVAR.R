library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)
source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/simu_data.R")
#source("./R/data_generator.R")

#sourceCpp("./src/cgVARSSL_dpe.cpp")
#sourceCpp("./src/cgVARSSL_dcpe.cpp")
sourceCpp("./src/cgSSL.cpp", rebuild = T)

q <- 10

graph <- g_model1(q)
#Omega <- graph$Omega
#Sigma <- graph$Sigma
#B <- graph$Omega
Sigma <- graph$Sigma
B <- diag(2*(c(1:q)%%2)-1)
#Sigma <- diag(q)

Y <- rCVARm(100,B,Sigma)

varSSL_dpe_res <- cgVARSSL_dpe(Y,lambdas = list(lambda1 = 1, lambda0 = seq(10, (nrow(Y)-1), length = 10)),
                    xis = list(xi1 = 0.01 * (nrow(Y)-1), xi0 = seq(0.1 * ((nrow(Y)-1)), (nrow(Y)-1), length = 10)),
                    theta_hyper_params = c(1, ncol(Y) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 10000,
                    eps = 1e-6,
                    s_max_condition = 10*(nrow(Y)-1),
                    obj_counter_max = 5,
                    verbose = 0)

varSSL_dcpe_res <- cgVARSSL_dcpe(Y,
                    lambdas = list(lambda1 = 1, lambda0 = seq(10, (nrow(Y)-1), length = 10)),
                    xis = list(xi1 = 0.01 * (nrow(Y)-1), xi0 = seq(0.1 * (nrow(Y)-1), (nrow(Y)-1), length = 10)),
                    theta_hyper_params = c(1, ncol(Y) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 10000,
                    eps = 1e-6,
                    #s_max_condition = 10*(nrow(Y)-1),
                    #obj_counter_max = 5,
                    verbose = 1)

error_B(varSSL_dpe_res$B,B)
error_Omega(varSSL_dpe_res$Omega, graph$Omega)

error_B(varSSL_dcpe_res$B,B)
error_Omega(varSSL_dcpe_res$Omega, graph$Omega)
