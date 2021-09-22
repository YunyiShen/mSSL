library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)



source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/simu_data.R")
sourceCpp("./src/fstarSSL.cpp",rebuild = TRUE)
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
mu <- 0*runif(q,-2,2)+2
mu <- 0 * mu 

Y <- rSTARlogfloor(X,B,mu,Sigma)
Y <- rSTARidfloor0(X,B,mu,Sigma)

lower <- logfloorlink(Y,FALSE)
upper <- logfloorlink(Y,TRUE)

star_res_dpe <- fstarmSSL(Y,X,link = logfloorlink,verbose = TRUE)
star_res_dpe <- fstarmSSL(Y,X,link = idfloor0link,verbose = TRUE)
