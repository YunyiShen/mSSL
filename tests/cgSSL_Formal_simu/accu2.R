library(CARlasso)
library(MASS)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("./src/cgSSL_dcpe.cpp")
sourceCpp("./src/cgSSL_dpe.cpp")
sourc_R <- lapply(list.files("./R",
                    pattern = "generator|error", 
                    full.names = TRUE), 
            source) # useful helpers

set.seed(12345)
file_num <- 2

res_loss_file <- "./tests/cgSSL_Formal_simu/Res/loss_"
res_graph_Omega_file <- "./tests/cgSSL_Formal_simu/Res/graph_Omega_"
res_graph_B_file <- "./tests/cgSSL_Formal_simu/Res/graph_B_"

res_loss_file <- paste0(res_loss_file,file_num,".csv")
res_graph_Omega_file <- paste0(res_graph_Omega_file,file_num,".csv")
res_graph_B_file <- paste0(res_graph_B_file,file_num,".csv")



#ns <- 100
ps <- c(20)
qs <- c(30)
ss <- c(.2)

nrep <- 50


sigma_models <- 3:4
retry <- 5

n_res_loss <- length(qs) * length(ss) * 
              length(ps) * length(sigma_models) * 
              nrep * 4 # 4 methods, dpe, dcpe, car, cara

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
    n <- 100
    for(p in ps){
        for(s in ss){

            B <- as.matrix( rsparsematrix(p, q, s, 
                        rand.x = function(n){runif(n,-2,2)}))
            X <- matrix(rnorm(n*p),n,p)
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
                    
                    cat("cgSSL-dpe\n")
                    for(jj in 1:retry){
                        res_dpe <- tryCatch( 
                            cgSSL_dpe(X,Y,
                                lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                                xis = list(xi1 = 0.01 * nrow(X), 
                                    xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)
                                ),
                                theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                                eta_hyper_params = c(1, ncol(Y)),
                                diag_penalty = 0,
                                max_iter = 10000,
                                eps = 1e-3,
                                s_max_condition = 10*nrow(X),
                                obj_counter_max = 5,
                                verbose = 0),
                            error = function(e){return(list())} )
                        if(length(res_dpe)>0) break
                        else cat("retry ",jj,"\n")
                    }
                    res_loss[i_res_loss,1:6] <- c(q,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "cgSSL-dpe"
                    res_loss$logL2B[i_res_loss] <- log_l2(B,res_dpe$B)
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


                    cat("cgSSL-dcpe\n")
                    for(jj in 1:retry){
                        res_dcpe <- tryCatch( 
                            cgSSL_dcpe(X,Y,
                                lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                                xis = list(xi1 = 0.01 * nrow(X), 
                                    xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)
                                ),
                                theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                                eta_hyper_params = c(1, ncol(Y)),
                                diag_penalty = 0,
                                max_iter = 10000,
                                eps = 1e-3,
                                verbose = 0),
                            error = function(e){return(list())} )
                        if(length(res_dcpe)>0) break
                        else cat("retry ",jj,"\n")
                    }
                    res_loss[i_res_loss,1:6] <- c(q,p,n,s,mod,i)
                    res_loss$algo[i_res_loss] <- "cgSSL-dcpe"
                    res_loss$logL2B[i_res_loss] <- log_l2(B,res_dcpe$B)
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

                    # start to use CARlasso 
                    df <- data.frame(Y, X)
                    colnames(df) <- c(paste0("y",1:q), paste0("x",1:p))
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
                    res_loss$logL2B[i_res_loss] <- log_l2(B,CAR_B)
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

                    # adaptive
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
                    res_loss$logL2B[i_res_loss] <- log_l2(B,CARA_B)
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

