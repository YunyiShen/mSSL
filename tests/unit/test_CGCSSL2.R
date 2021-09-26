library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)
source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/plot_support.R")
#sourceCpp("./src/cgSSL_dcpe.cpp")
#sourceCpp("./src/cgSSL_dpe.cpp")
sourceCpp("./src/cgSSL.cpp",rebuild = TRUE)

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
    Y[i,] <- mvrnorm(n = 1, Sigma %*% XB[i,], Sigma)
}


GCSSL_dpe_res <- cgSSL_dpe(X,Y,lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                    xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                    theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 10000,
                    eps = 1e-6,
                    s_max_condition = 10*nrow(X),
                    obj_counter_max = 5,
                    verbose = 1)

GCSSL_dcpe_res <- cgSSL_dcpe(X,Y,
                    lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                    xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                    theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                    eta_hyper_params = c(1, ncol(Y)),
                    diag_penalty = 0,
                    max_iter = 10000,
                    eps = 1e-3,
                    #s_max_condition = 10*nrow(X),
                    #obj_counter_max = 5,
                    verbose = 1)

resp <- paste0(paste0("y",1:q), collapse = "+")
pred <- paste0(paste0("x",1:p), collapse = "+")
formu <- as.formula( paste0(resp, "~", pred) )

df <- data.frame(Y, X)
colnames(df) <- c(paste0("y",1:q), paste0("x",1:p))

baseline <- CARlasso::CARlasso(formu, data = df, adaptive = TRUE, n_iter = 5000)
baseline <- CARlasso::horseshoe(baseline)

CARlasso:::stein_loss(graph$Omega, GCSSL_dcpe_res$Omega)
CARlasso:::stein_loss(graph$Omega, GCSSL_dpe_res$Omega)
CARlasso:::stein_loss(graph$Omega, baseline$point_est$Omega)


pdf("./one_exp.pdf", width = 8, height = 5)
par(mfrow = c(2,4))
image(B, main = "ground truth", ylab = "B")
image(GCSSL_dpe_res$B, main = "cgSSL DPE")
image(GCSSL_dcpe_res$B, main = "cgSSL DCPE")
image(baseline$point_est$beta, main = "CAR-ALASSO")

image(Omega, ylab = "Omega")
image(GCSSL_dpe_res$Omega)
image(GCSSL_dcpe_res$Omega)
image(baseline$point_est$Omega)
dev.off()

erB_dcpe <- error_B(GCSSL_dcpe_res$B, B)
erB_dpe <- error_B(GCSSL_dpe_res$B, B)
erB_carlasso <- error_B(baseline$horseshoe_binary$B_binary * baseline$point_est$beta, B)

erB <- rbind(erB_dcpe,erB_dpe, erB_carlasso)
rownames(erB) <- c("cgSSL_dcpe","cgSSL_dpe","CAR_Alasso")

erO_dcpe <- error_Omega(GCSSL_dcpe_res$Omega, Omega)
erO_dpe <- error_Omega(GCSSL_dpe_res$Omega, Omega)
erO_carlasso <- error_Omega(baseline$horseshoe_binary$Omega_binary * baseline$point_est$Omega, 
    Omega)

erO <- rbind(erO_dcpe,erO_dpe, erO_carlasso)
rownames(erO) <- c("cgSSL_dcpe","cgSSL_dpe","CAR_Alasso")

erB
erO

plot_support(GCSSL_dpe_res)
