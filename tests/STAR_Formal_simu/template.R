
args <- commandArgs(trailingOnly=TRUE)
library(MASS)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(mSSL)

rSTARlogfloor <- function(X,B,mu,Sigma){
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    res <- XB+E
    res <- floor(exp(res))
    return(res)
}

rSTARidfloor0 <- function(X,B,mu,Sigma){
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    res <- XB+E
    res <- floor(res)
    return(res * (res>=0))
}

rSTARsqrtfloor0 <- function(X,B,mu,Sigma){
    XB <- X %*% B
    E <- mvrnorm(nrow(X),mu,Sigma)
    res <- XB+E
    temp <- floor(res^2)
    return(temp * (res>=0))
}


g_model1 <- function(k, rho=.7, tol = 1e-10){
  temp <- matrix(rep(1:k,k),ncol = k)
  Sigma <- rho ^ (abs(temp-t(temp)))
  Omega <- solve(Sigma)
  Omega <- Omega * (abs(Omega)>tol)
  return(list(Sigma = Sigma, Omega = Omega))
}







file_num <- paste0(args[1],"-",args[2])
base_dir <- paste0("./res_",args[1],"-",args[2])

res_loss_file <- paste0(base_dir,"/loss_")
res_graph_Omega_file <- paste0(base_dir, "/graph_Omega_")
res_graph_B_file <- paste0(base_dir, "/graph_B_")

res_loss_file <- paste0(res_loss_file,file_num,".csv")
res_graph_Omega_file <- paste0(res_graph_Omega_file,file_num,".csv")
res_graph_B_file <- paste0(res_graph_B_file,file_num,".csv")


p <- 10
q <- 5
s <- .3
n <- 100
n_test <- 25


all_transformations <- c("id0","log","sqrt")
sampling_functions <- list(rSTARidfloor0, rSTARlogfloor, rSTARsqrtfloor0)


all_methods <- c("id0","log","sqrt","cdf-all","cdf-ind","mSSL")

links_first_4 <- list(idfloor0link, logfloorlink, sqrtfloor0link, KowalWufloorlink)

internal_procss_id <- floor( as.numeric( args[2] )/18) # this determined which repeat
internal_experiment_id <- as.numeric( args[2] ) %% 18

transformation_truth <- (internal_experiment_id %% 3) + 1 # which data generating process to do 
method_used <- floor( internal_experiment_id/3 ) + 1
retry <- 5

n_res_loss <- 1

res_loss <- data.frame(matrix(NA,nrow = n_res_loss,ncol = 8))
colnames(res_loss) <- c("q","p","n","s","mod",
                        "irep","algo","logMSE")

res_graph_Omega <- data.frame(matrix(NA,nrow = n_res_loss,ncol = 18))
colnames(res_graph_Omega) <- c("q","p","n","s","mod",
                               "irep","algo","TP", "TN", 
                               "FP", "FN", "SEN", "SPE", 
                               "PREC", "ACC",  "F1", 
                               "MCC", "FROB")

res_graph_B <- data.frame(matrix(NA,nrow = n_res_loss,ncol = 18))
colnames(res_graph_B) <- c("q","p","n","s","mod",
                           "irep","algo","TP", "TN", 
                           "FP", "FN", "SEN", "SPE", 
                           "PREC", "ACC",  "F1", 
                           "MCC", "FROB")


set.seed(42)
B <- as.matrix( rsparsematrix(p, q, s, 
                                    rand.x = function(n){runif(n,-2,2)}))
X <- matrix(rnorm((n + n_test) * p),ncol = p)

mu <- rep(2,q)

set.seed(as.numeric(args[1])+internal_procss_id )
graph_tmp <- g_model1(k=q)
Omega <- graph_tmp$Omega
Sigma <- graph_tmp$Sigma

## sampling code goes here
Y <- sampling_functions[[transformation_truth]](X = X, B=B, Sigma = Sigma, mu = mu)
##

X_test <- X[-c(1:n),]
X <- X[1:n,]

Y_test <- Y[-c(1:n),]
Y <- Y[1:n,]
          


if(method_used %in% c(1,2,3,4) ){ 
  
  for(jj in 1:retry){
      res <- tryCatch( 
          fstarmSSL(X=X, Y=Y, link = links_first_4[[method_used]], verbose = FALSE),
            error = function(e){return(list())} )
          if(length(res)>0) break
          else cat("retry ",jj,"\n")
  }

  pred_mean <- X_test %*% res$B  + rep(1,nrow(X_test))%*% t(res$alpha)
  pred_mean_trans <- matrix( res$link_res$invtransforming(pred_mean), nrow = nrow(Y_test) )        
}

if(method_used==5){
  for(jj in 1:retry){
      res <- tryCatch( 
          fstarmSSL(X=X, Y=Y, link = KowalWufloorlink, verbose = FALSE, overall = FALSE),
            error = function(e){return(list())} )
          if(length(res)>0) break
          else cat("retry ",jj,"\n")
  }

  pred_mean <- X_test %*% res$B  + rep(1,nrow(X_test))%*% t(res$alpha)
  pred_mean_trans <- lapply(1:q, function(i, pred_mean, transforms){
    transforms[[i]](pred_mean[,i])
  }, pred_mean, res$link_res$invtransforming)
  pred_mean_trans <- Reduce(cbind, pred_mean_trans)
}

if(method_used==6){
    for(jj in 1:retry){
      res <- tryCatch( 
          mSSL(X=X, Y=Y, verbose = FALSE, cg = FALSE),
            error = function(e){return(list())} )
          if(length(res)>0) break
          else cat("retry ",jj,"\n")
  }
  pred_mean_trans <- X_test %*% res$B  + rep(1,nrow(X_test))%*% t(res$alpha)
}

res_loss[1,1:6] <- c(q,p,n,s,transformation_truth,method_used)
res_loss$algo[1] <- all_methods[method_used]
res_loss$logMSE[1] <- log(sum((pred_mean_trans-Y_test)^2)/n_test)  
res_loss$mod[1] <- all_transformations[res_loss$mod[1]]
write.csv(res_loss,res_loss_file)

          
res_graph_Omega[1,1:6] <- c(q,p,n,s,transformation_truth,internal_experiment_id)
res_graph_B[1,1:6] <- c(q,p,n,s,transformation_truth,internal_experiment_id)
res_graph_Omega[1,7] <- all_methods[method_used]
res_graph_B[1,7] <- all_methods[method_used]
res_graph_Omega[1,8:18] <- error_Omega(res$Omega, Omega)
res_graph_B[1,8:18] <- error_B(res$B, B)

res_graph_Omega$mod[1] <- all_transformations[res_graph_Omega$mod[1]]
res_graph_B$mod[1] <- all_transformations[res_graph_B$mod[1]]


write.csv(res_graph_B,res_graph_B_file)
write.csv(res_graph_Omega,res_graph_Omega_file)

