library(Rcpp)
library(RcppArmadillo)
library(MASS)
library(Matrix)

KowalWufloorlink2 <- function(Y, distribution = "np", grid = NULL){
    p <- ncol(Y)
    trans_funs <- g_cdf(Y, distribution)

    lower <- lapply(1:p, function(i,Y,trans_funs){
        trans_funs(Y[,i])
    },Y,trans_funs)
    lower <- Reduce(cbind, lower)
    upper <- lapply(1:p, function(i,Y,trans_funs){
        trans_funs(Y[,i] + 1)
    },Y,trans_funs)
    upper <- Reduce(cbind, upper)
    if(is.null(grid)){
        grid <- seq(min(5*min(lower),-1),max(5*max(upper),1),by = .1)
    }
    invtransforming <-  g_inv_approx(trans_funs, grid)

    return(list(lower = lower, upper = upper, invtransforming = invtransforming))

}


source("./R/graph_generator.R")
source("./R/error_B.R")
source("./R/error_Omega.R")
source("./R/simu_data.R")
sourceCpp("./src/mSSL.cpp", rebuild = T)
#sourceCpp("./src/fstarSSL_dpe.cpp",rebuild = TRUE)
#sourceCpp("./src/fstarSSL_dcpe.cpp")
source("./R/starSSL.R")

set.seed(12345)
p <- 3
q <- 5
n <- 300
B <- as.matrix( rsparsematrix(p, q, .3, rand.x = function(n){runif(n,-2,2)}))

X <- matrix(rnorm(p * n), nrow = n, ncol = p)

graph <- g_model1(q)
Sigma <- graph$Sigma
Omega <- graph$Omega
mu <- 0*runif(q,-2,2)+2
mu <- 0 * mu + 2

Y <- rSTARlogfloor(X,B,mu,Sigma)
Y <- rSTARidfloor(X,B,mu,Sigma)
Y <- rSTAR(X,B,mu,Sigma, link = function(x){
    temp <- floor(x^2)
    temp[x<=0] <- 0
    return(temp)
})

bounds <- sqrtfloor0link(Y)
bounds2 <- KowalWufloorlink(Y)
par(mfrow = c(2,3))

for(i in 1:5){
t_lower <- bounds2$lower[,i] 
t_upper <- bounds2$upper[,i] 
plot(Y[,i],bounds2$invtransforming[[i]](t_upper), type = "p", col = "red")  
points(Y[,i],bounds2$invtransforming[[i]](t_lower), col = "blue")
lines(Y[,i],Y[,i])   
} 

par(mfrow = c(2,3))
for(i in 1:5){
    xs <- seq(-10,10,0.1)
    sqrx <- xs^2
    sqrx[xs<=0] <- 0

    invtrsres <- bounds2$invtransforming[[i]](xs)
    plot(xs, invtrsres, type = "l", col = "red")
    lines(xs, sqrx, col = "blue")
    lines(xs, g_cdf(Y[,i])(xs))
}


star_res_dpe <- fstarmSSL(Y,X,link = KowalWufloorlink,verbose = TRUE, eps=1e-1)
star_res_dpe1 <- fstarmSSL(Y,X,link = sqrtfloor0link,verbose = TRUE, eps=1e-1, s_max_condition = 1e4)

pred_mean <- X %*% star_res_dpe$B  + rep(1,nrow(X))%*% t(star_res_dpe$alpha)
pred_mean_trans <- matrix( star_res_dpe$link_res$invtransforming(pred_mean), nrow = nrow(Y) )

pred_mean1 <- X %*% star_res_dpe1$B  + rep(1,nrow(X))%*% t(star_res_dpe1$alpha)
pred_mean_trans1 <- matrix( star_res_dpe1$link_res$invtransforming(pred_mean1), nrow = nrow(Y) )

mSSL_res_dpe <- mSSL::mSSL(Y = Y, X = X)
pred_mean_mssl <- X %*% mSSL_res_dpe$B  + rep(1,nrow(X))%*% t(mSSL_res_dpe$alpha)


