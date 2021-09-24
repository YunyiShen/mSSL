library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)



source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/simu_data.R")
sourceCpp("./src/cgfstarSSL_dpe.cpp",rebuild = TRUE)
sourceCpp("./src/cgfstarSSL_dcpe.cpp",rebuild = TRUE)
source("./R/starSSL.R")


set.seed(12345)
p <- 2
q <- 5
n <- 50
B <- as.matrix( rsparsematrix(p, q, 0.3, rand.x = function(n){runif(n,-2,2)}))

X <- matrix(rnorm(p * n), nrow = n, ncol = p)

graph <- g_model1(q)
Sigma <- graph$Sigma
Omega <- graph$Omega
mu <- 0*runif(q,-2,2)+1


Y <- rcgSTAR(X,B,mu,Sigma)

lower <- logfloorlink(Y,FALSE)
upper <- logfloorlink(Y,TRUE)

cgstar_res_dpe <- fstarmSSL(Y,X,link = logfloorlink,cg = TRUE,verbose = TRUE)
cgstar_res_dcpe <- fstarmSSL(Y,X,link = logfloorlink, cg = TRUE, condexp = TRUE,verbose = TRUE, eps = 1e-3)

