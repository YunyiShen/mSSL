library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)



source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/simu_data.R")
sourceCpp("./src/mSSL.cpp", rebuild = T)
sourceCpp("./src/fstarSSL_dpe.cpp",rebuild = TRUE)
sourceCpp("./src/fstarSSL_dcpe.cpp")
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
Y <- rSTARidfloor(X,B,mu,Sigma)

bounds <- logfloorlink(Y)
bounds2 <- KowalWufloorlink(Y)
t_lower <- bounds2$lower[,1] 
t_upper <- bounds2$upper[,1] 
plot(Y[,1],bounds2$invtransforming[[1]](t_upper))  
points(Y[,1],bounds2$invtransforming[[1]](t_lower))
lines(Y,Y)    


star_res_dpe <- fstarmSSL(Y,X,link = idfloorlink,verbose = TRUE, eps=1e-2)

star_res_dpe1 <- fstarmSSL(Y,X,link = KowalWufloorlink,verbose = TRUE, eps=1e-2, s_max_condition = 1e4)
