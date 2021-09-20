library(Rcpp)
library(RcppArmadillo)
library(MASS)

sourceCpp("./src/cgquic.cpp")
sourceCpp("./src/sandbox.cpp")
sourceCpp("./src/quic_ori.cpp")
source("./R/graph_generator.R")

q <- 10
n <- 500

set.seed(
    12345
)
mus <- matrix(rnorm(q*n),n,q)
M <- t(mus) %*% mus / n

graph <- g_model1(q) # get Sigma and Omega
Omega <- graph$Omega
Sigma <- graph$Sigma

Y <- mus # just for size match

for(i in 1:n){
    Y[i,] <- MASS::mvrnorm(n=1, mu =  mus[i,] %*% Sigma, Sigma)
}


S <- t(Y) %*% Y / n


D <- diag(q)
D <- S
Rho <- matrix(10,q,q)/n
diag(Rho) <- 0
Rho1 <- Rho
#Rho1[upper.tri(Rho)] <- 0
#diag(Rho) <- 0
U <- D %*% Sigma
Q <- M %*% Sigma
normD <- 0.0
diffD <- 0.0

CoordinateDescentUpdate_ori(q, S, Rho, Omega, Sigma, U, D, 0,3, normD, diffD)
CoordinateDescentUpdate(q, S, Rho, Omega, Sigma, U, Q, D, 0,3, normD, diffD)
tol <- 1e-3

 
tt<-quic_res <- my_cgquic(q,S,M,Rho, tol, 2000)
tt2<-quic_res <- my_quic(q = q,S = cov(S),Rho = Rho, tol = tol, max_iter_quic = 2000)
CARlasso:::stein_loss(Omega, tt[,,1])  
CARlasso:::stein_loss(Omega, tt2[,,1])  

log(det(quic_res[,,1]))-sum(diag(S %*% quic_res[,,1])) - sum(diag(M %*% quic_res[,,2])) - sum(abs(Rho1 * quic_res[,,1])) 
log(det(graph$Omega))-sum(diag(S %*% graph$Omega)) - sum(diag(M %*% graph$Sigma)) - sum(abs(Rho1 * graph$Omega)) 
a = 0.01
log(det(quic_res[,,1]+a * graph$Omega))-sum(diag(S %*% (quic_res[,,1]+a * graph$Omega))) - sum(diag(M %*% solve(quic_res[,,1]+a * graph$Omega))) - sum(abs(Rho1 * (quic_res[,,1]+a*graph$Omega))) 
