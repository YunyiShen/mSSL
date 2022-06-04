g_model1 <- function(k, rho=.7, tol = 1e-10){
  temp <- matrix(rep(1:k,k),ncol = k)
  Sigma <- rho ^ (abs(temp-t(temp)))
  Omega <- solve(Sigma)
  Omega <- Omega * (abs(Omega)>tol)
  return(list(Sigma = Sigma, Omega = Omega))
}


g_model2 <- function(k, rho = .5){
  temp <- matrix(rep(1:k,k),ncol = k)
  Omega <- rho ^ (abs(temp-t(temp))) * (abs(temp-t(temp)) <= 2)
  return(list(Omega = Omega, Sigma = solve(Omega)))
}

g_model3 <- function(k, rho = .5){
  row_ind <- matrix(rep(1:k,k),ncol = k)
  col_ind <- t(row_ind)
  
  Sigma <- diag(1,k,k)
  
  Sigma[(1<=row_ind ) &
          (1<=col_ind ) &
          (row_ind != col_ind) & 
          (row_ind <= k/2) & 
          (col_ind <= k/2)] <- rho
  Sigma[((1+k/2)<=row_ind ) & 
          ((1+k/2)<=col_ind ) &
          (row_ind != col_ind)] <- rho # & 
  #(row_ind <= 10) & 
  #(col_ind <= 10)] 
  Omega <- solve(Sigma)
  
  return(list(Sigma = Sigma, Omega = Omega))
}

g_model4 <- function(k, rho=.1){
  
  Omega <- diag(1,k,k)
  Omega[2:k,1] <- rho
  Omega[1,2:k] <- rho
  Sigma <- solve(Omega)
  return(list(Sigma = Sigma, Omega = Omega))
}

g_model5 <- function(k, rhos = c(2,1,.9)){
  row_ind <- matrix(rep(1:k,k),ncol = k)
  col_ind <- t(row_ind)
  Omega <- diag(rhos[1],k,k)
  Omega[abs(row_ind-col_ind)==1] <- rhos[2]
  Omega[1,k] <- Omega[k,1] <- rhos[3]
  return(list(Sigma = solve(Omega), Omega = Omega))
}

g_model6 <- function(k, rhos = c(2,1)){
  Omega <- matrix(rhos[2], k, k)
  diag(Omega) <- rhos[1]
  
  return(list(Sigma = solve(Omega),Omega = Omega))
}


args <- commandArgs(trailingOnly=TRUE)
library(CARlasso)
library(MASS)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(mSSL)



file_num <- paste0(args[1],"-",args[2])
base_dir <- paste0("./res_",args[1],"-",args[2])

res_loss_file <- paste0(base_dir,"/loss_")
res_graph_Omega_file <- paste0(base_dir, "/graph_Omega_")
res_graph_B_file <- paste0(base_dir, "/graph_B_")

res_loss_file <- paste0(res_loss_file,file_num,".csv")
res_graph_Omega_file <- paste0(res_graph_Omega_file,file_num,".csv")
res_graph_B_file <- paste0(res_graph_B_file,file_num,".csv")

ps <- c(100)
qs <- c(30)
ss <- c(.2)

nrep <- 1


sigma_models <- c(as.numeric(args[3]))
retry <- 5

n_res_loss <- length(qs) * length(ss) * 
  length(ps) * length(sigma_models) * 
  nrep # 6 methods, dpe, dcpe, car, cara, with mSSL-dpe and mSSL-dcpe

res_loss <- data.frame(matrix(NA,nrow = n_res_loss,ncol = 9))
colnames(res_loss) <- c("q","p","n","s","mod",
                        "irep","algo","logL2B",
                        "steinOmega")

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
i_res_loss <- 1
i_res_graph_Omega <- 1
i_res_graph_B <- 1
for(q in qs) {
  n <- 400
  for(p in ps){
    for(s in ss){
      set.seed(42)
      B <- as.matrix( rsparsematrix(p, q, s, 
                                    rand.x = function(n){runif(n,-2,2)}))
      X <- matrix(rnorm(n*p),n,p)
      set.seed(as.numeric(args[1])+floor(as.numeric(args[2])/6))
	resp <- paste0(paste0("y",1:q), collapse = "+")
      pred <- paste0(paste0("x",1:p), collapse = "+")
      formu <- as.formula( paste0(resp, "~", pred) )
      
      for(mod in sigma_models){
        graph_tmp <- do.call(paste0("g_model",mod),list(k=q))
        Omega <- graph_tmp$Omega
        Sigma <- graph_tmp$Sigma
        
	for(i in 1:nrep){
          cat("q =",q,",p =",p,",s =",s,",mod =",mod,",rep =",i,"\n")
          Y <- rCARm(n,B,X,Sigma)
          if(as.numeric(args[2])%%6==0){
          cat("cgSSL-dpe\n")
          for(jj in 1:retry){
            res_dpe <- tryCatch( 
              mSSL(X=X,Y=Y,cg = TRUE,condexp = FALSE, verbose = FALSE),
              error = function(e){return(list())} )
            if(length(res_dpe)>0) break
            else cat("retry ",jj,"\n")
          }
          res_loss[i_res_loss,1:6] <- c(q,p,n,s,mod,i)
          res_loss$algo[i_res_loss] <- "cgSSL-dpe"
          res_loss$logL2B[i_res_loss] <- mSSL:::log_l2(B,res_dpe$B)
          res_loss$steinOmega[i_res_loss] <- CARlasso:::stein_loss_cpp(res_dpe$Omega, Omega)
          write.csv(res_loss,res_loss_file)
          i_res_loss <- i_res_loss + 1
          
          res_graph_Omega[i_res_graph_Omega,1:6] <- c(q,p,n,s,mod,i)
          res_graph_B[i_res_graph_B,1:6] <- c(q,p,n,s,mod,i)
          res_graph_B[i_res_graph_B,7] <- "cgSSL-dpe"
          res_graph_Omega[i_res_graph_Omega,7] <- "cgSSL-dpe"
          
          res_graph_Omega[i_res_graph_Omega,8:18] <- error_Omega(res_dpe$Omega, Omega)
          res_graph_B[i_res_graph_B,8:18] <- error_B(res_dpe$B, B)
          i_res_graph_Omega <- i_res_graph_Omega + 1
          i_res_graph_B <- i_res_graph_B + 1
          write.csv(res_graph_B,res_graph_B_file)
          write.csv(res_graph_Omega,res_graph_Omega_file)
          
          }
          
	  if(as.numeric(args[2])%%6==1){
	  cat("cgSSL-dcpe\n")
          for(jj in 1:retry){
            res_dcpe <- tryCatch( 
              mSSL(X=X,Y=Y,cg = TRUE,condexp = TRUE, verbose = FALSE),
              error = function(e){return(list())} )
            if(length(res_dcpe)>0) break
            else cat("retry ",jj,"\n")
          }
          res_loss[i_res_loss,1:6] <- c(q,p,n,s,mod,i)
          res_loss$algo[i_res_loss] <- "cgSSL-dcpe"
          res_loss$logL2B[i_res_loss] <- mSSL:::log_l2(B,res_dcpe$B)
          res_loss$steinOmega[i_res_loss] <- CARlasso:::stein_loss_cpp(res_dcpe$Omega, Omega)
          write.csv(res_loss,res_loss_file)
          i_res_loss <- i_res_loss + 1
          
          res_graph_Omega[i_res_graph_Omega,1:6] <- c(q,p,n,s,mod,i)
          res_graph_B[i_res_graph_B,1:6] <- c(q,p,n,s,mod,i)
          res_graph_Omega[i_res_graph_Omega,7] <- "cgSSL-dcpe"
          res_graph_B[i_res_graph_B,7] <- "cgSSL-dcpe"
          res_graph_Omega[i_res_graph_Omega,8:18] <- error_Omega(res_dcpe$Omega, Omega)
          res_graph_B[i_res_graph_B,8:18] <- error_B(res_dcpe$B, B)
          i_res_graph_Omega <- i_res_graph_Omega + 1
          i_res_graph_B <- i_res_graph_B + 1
          write.csv(res_graph_B,res_graph_B_file)
          write.csv(res_graph_Omega,res_graph_Omega_file)
          
          # use of mSSL
          }
	  if(as.numeric(args[2])%%6==2){
          cat("mSSL-dpe\n")
          for(jj in 1:retry){
            res_mdpe <- tryCatch( 
              mSSL(X=X,Y=Y,cg = FALSE,verbose = FALSE ),
              error = function(e){return(list())} )
            if(length(res_mdpe)>0) break
            else cat("retry ",jj,"\n")
          }
          res_loss[i_res_loss,1:6] <- c(q,p,n,s,mod,i)
          res_loss$algo[i_res_loss] <- "mSSL-dpe"
          res_loss$logL2B[i_res_loss] <- mSSL:::log_l2(B,res_mdpe$B %*% res_mdpe$Omega)
          res_loss$steinOmega[i_res_loss] <- CARlasso:::stein_loss_cpp(res_mdpe$Omega, Omega)
          write.csv(res_loss,res_loss_file)
          i_res_loss <- i_res_loss + 1
          
          res_graph_Omega[i_res_graph_Omega,1:6] <- c(q,p,n,s,mod,i)
          res_graph_B[i_res_graph_B,1:6] <- c(q,p,n,s,mod,i)
          res_graph_B[i_res_graph_B,7] <- "mSSL-dpe"
          res_graph_Omega[i_res_graph_Omega,7] <- "mSSL-dpe"
          
          res_graph_Omega[i_res_graph_Omega,8:18] <- error_Omega(res_mdpe$Omega, Omega)
          res_graph_B[i_res_graph_B,8:18] <- error_B(res_mdpe$B%*%res_mdpe$Omega, B)
          i_res_graph_Omega <- i_res_graph_Omega + 1
          i_res_graph_B <- i_res_graph_B + 1
          write.csv(res_graph_B,res_graph_B_file)
          write.csv(res_graph_Omega,res_graph_Omega_file)
          
          }
	  if(as.numeric(args[2])%%6==3){
          cat("mSSL-dcpe\n")
          for(jj in 1:retry){
            res_mdcpe <- tryCatch( 
              mSSL(X=X,Y=Y,cg = FALSE,condexp = TRUE, verbose = FALSE),
              error = function(e){return(list())} )
            if(length(res_mdcpe)>0) break
            else cat("retry ",jj,"\n")
          }
          res_loss[i_res_loss,1:6] <- c(q,p,n,s,mod,i)
          res_loss$algo[i_res_loss] <- "mSSL-dcpe"
          res_loss$logL2B[i_res_loss] <- mSSL:::log_l2(B,res_mdcpe$B%*%res_mdcpe$Omega)
          res_loss$steinOmega[i_res_loss] <- CARlasso:::stein_loss_cpp(res_mdcpe$Omega, Omega)
          write.csv(res_loss,res_loss_file)
          i_res_loss <- i_res_loss + 1
          
          res_graph_Omega[i_res_graph_Omega,1:6] <- c(q,p,n,s,mod,i)
          res_graph_B[i_res_graph_B,1:6] <- c(q,p,n,s,mod,i)
          res_graph_Omega[i_res_graph_Omega,7] <- "mSSL-dcpe"
          res_graph_B[i_res_graph_B,7] <- "mSSL-dcpe"
          res_graph_Omega[i_res_graph_Omega,8:18] <- error_Omega(res_mdcpe$Omega, Omega)
          res_graph_B[i_res_graph_B,8:18] <- error_B(res_mdcpe$B%*%res_mdcpe$Omega, B)
          i_res_graph_Omega <- i_res_graph_Omega + 1
          i_res_graph_B <- i_res_graph_B + 1
          write.csv(res_graph_B,res_graph_B_file)
          write.csv(res_graph_Omega,res_graph_Omega_file)
          }
          # start to use CARlasso 
          df <- data.frame(Y, X)
          colnames(df) <- c(paste0("y",1:q), paste0("x",1:p))
          if(as.numeric(args[2])%%6==4){
	  cat("CARlasso\n")
          for(jj in 1:retry){
            res_CAR <- tryCatch( 
              CARlasso(formu, data = df, n_iter = 5000, 
                       thin_by = 5,verbos = FALSE),
              error = function(e){return(list())} )
            if(length(res_CAR)>0) break
            else cat("retry ",jj,"\n")
          }
          res_CAR <- CARlasso::horseshoe(res_CAR)
          CAR_B <- res_CAR$point_est$beta * res_CAR$horseshoe_binary$B_binary
          CAR_Omega <- res_CAR$point_est$Omega * res_CAR$horseshoe_binary$Omega_binary
          
          res_loss[i_res_loss,1:6] <- c(q,p,n,s,mod,i)
          res_loss$algo[i_res_loss] <- "CAR"
          res_loss$logL2B[i_res_loss] <- mSSL:::log_l2(B,CAR_B)
          res_loss$steinOmega[i_res_loss] <- CARlasso:::stein_loss_cpp(CAR_Omega, Omega)
          write.csv(res_loss,res_loss_file)
          i_res_loss <- i_res_loss + 1
          
          res_graph_Omega[i_res_graph_Omega,1:6] <- c(q,p,n,s,mod,i)
          res_graph_B[i_res_graph_B,1:6] <- c(q,p,n,s,mod,i)
          res_graph_Omega[i_res_graph_Omega,7] <- "CAR"
          res_graph_B[i_res_graph_B,7] <- "CAR"
          res_graph_Omega[i_res_graph_Omega,8:18] <- error_Omega(CAR_Omega, Omega)
          res_graph_B[i_res_graph_B,8:18] <- error_B(CAR_B, B)
          i_res_graph_Omega <- i_res_graph_Omega + 1
          i_res_graph_B <- i_res_graph_B + 1
          write.csv(res_graph_B,res_graph_B_file)
          write.csv(res_graph_Omega,res_graph_Omega_file)
          }
          # adaptive
          if(as.numeric(args[2])%%6==5){
	  cat("CARlasso-A\n")
          for(jj in 1:retry){
            res_CARA <- tryCatch( 
              CARlasso(formu, data = df, n_iter = 5000, 
                       thin_by = 5,verbos = FALSE, adaptive = TRUE),
              error = function(e){return(list())} )
            if(length(res_CARA)>0) break
            else cat("retry ",jj,"\n")
          }
          res_CARA <- CARlasso::horseshoe(res_CARA)
          CARA_B <- res_CARA$point_est$beta * res_CARA$horseshoe_binary$B_binary
          CARA_Omega <- res_CARA$point_est$Omega * res_CARA$horseshoe_binary$Omega_binary
          
          res_loss[i_res_loss,1:6] <- c(q,p,n,s,mod,i)
          res_loss$algo[i_res_loss] <- "CAR-A"
          res_loss$logL2B[i_res_loss] <- mSSL:::log_l2(B,CARA_B)
          res_loss$steinOmega[i_res_loss] <- CARlasso:::stein_loss_cpp(CARA_Omega, Omega)
          write.csv(res_loss,res_loss_file)
          i_res_loss <- i_res_loss + 1
          
          res_graph_Omega[i_res_graph_Omega,1:6] <- c(q,p,n,s,mod,i)
          res_graph_B[i_res_graph_B,1:6] <- c(q,p,n,s,mod,i)
          res_graph_Omega[i_res_graph_Omega,7] <- "CAR-A"
          res_graph_B[i_res_graph_B,7] <- "CAR-A"
          res_graph_Omega[i_res_graph_Omega,8:18] <- error_Omega(CARA_Omega, Omega)
          res_graph_B[i_res_graph_B,8:18] <- error_B(CARA_B, B)
          i_res_graph_Omega <- i_res_graph_Omega + 1
          i_res_graph_B <- i_res_graph_B + 1
          write.csv(res_graph_B,res_graph_B_file)
          write.csv(res_graph_Omega,res_graph_Omega_file)
        }
	}
      }
    }
  }
}

