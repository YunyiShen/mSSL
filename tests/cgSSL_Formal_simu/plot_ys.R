get_result_list <- function(Psi_results, Omega_results){
  res <- list()
  algo <- c("LASSO", "CAR", "CAR-A", "cgSSL-dcpe", "cgSSL-dpe")
  outcomes <- c("SEN", "PREC","FROB")
  for(out in outcomes){
    tmp <- data.frame("LASSO" = rep(NA, times = 501),
                      "CAR" = rep(NA, times = 501),
                      "CAR-A" = rep(NA, times = 501),
                      "cgSSL-dcpe" = rep(NA, times = 501),
                      "cgSSL-dpe" = rep(NA, times = 501))
    for(i in 1:5){
      index <- which(Psi_results[,"algo"] == algo[i])
      #if(algo[i] != "LASSO") tmp[,i] <- Psi_results[index,out]
      #else tmp[1:500,i] <- Psi_results[index, out]
      tmp[1:length(index),i] <- Psi_results[index, out]
    }
    `[[`(res,paste0(out, "_psi_results")) <- tmp
    
    tmp <- data.frame("LASSO" = rep(NA, times = 501),
                      "CAR" = rep(NA, times = 501),
                      "CAR-A" = rep(NA, times = 501),
                      "cgSSL-dcpe" = rep(NA, times = 501),
                      "cgSSL-dpe" = rep(NA, times = 501))
    for(i in 1:5){
      index <- which(Omega_results[,"algo"] == algo[i])
      #if(algo[i] != "LASSO") tmp[,i] <- Omega_results[index,out]
      #else tmp[1:500,i] <- Omega_results[index, out]
      tmp[1:length(index),i] <- Omega_results[index,out]
    }
    `[[`(res,paste0(out, "_omega_results")) <- tmp
    
  }
  return(res)
}

plot_result_list <- function(result_list, modelname, title=TRUE){
  boxplot(result_list$SEN_psi_results, horizontal = TRUE, main = ifelse(title, expression(Psi~"sensitivity"), ""),
          medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
          ylim = c(0,1), yaxt = "n", ylab = modelname)
  text(x = 0.2, y = 1, labels = "cgLASSO", cex = 1.2)
  text(x = 0.2, y = 2, labels = "CAR", cex = 1.2)
  text(x = 0.2, y = 3, labels = "CAR-A", cex = 1.2)
  text(x = 0.22, y = 4, labels = "cgSSL-DCPE", cex = 1.2)
  text(x = 0.2, y = 5, labels = "cgSSL-DPE", cex = 1.2)
  boxplot(result_list$PREC_psi_results, horizontal = TRUE, main = ifelse(title, expression(Psi~"precision"),""),
          medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
          ylim = c(0,1), yaxt = "n")
  boxplot(result_list$FROB_psi_results, horizontal = TRUE, main = ifelse(title, expression(Psi~"Frobenius"),""),
          medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
           yaxt = "n")
  
  boxplot(result_list$SEN_omega_results, horizontal = TRUE, main = ifelse(title, expression(Omega~"sensitivity"),""),
          medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
          ylim = c(0,1), yaxt = "n")
  
  boxplot(result_list$PREC_omega_results, horizontal = TRUE, main = ifelse(title, expression(Omega~"precision"),""),
          medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
          ylim = c(0,1), yaxt = "n")
  boxplot(result_list$FROB_omega_results, horizontal = TRUE, main = ifelse(title, expression(Omega~"Frobenius"),""),
          medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
           yaxt = "n")
}


get_result_table <- function(res_list){
  mean_sd <- lapply(res_list, function(w){
    as.matrix(w) |>
     apply(2,function(ww){
        paste0(round(mean(ww,na.rm = T),2)," (",round(sd(ww,na.rm = T),2),")")
      }) 
  }) |> Reduce(f = cbind)
  colnames(mean_sd) <- NULL
  #rownames(mean_sd) <- paste0("\\textit{",rownames(mean_sd),"}")
  cat(knitr::kable(mean_sd[,c(1,3,5,2,4,6)], "latex",booktabs = TRUE),"\n", escape = F)
}

p <- 20
q <- 30
n <- 100

prefix <- paste0("p", p, "q", q, "n", n)
Psi_file <- file.path(prefix, paste0(prefix, "_graph_B_full.csv"))
Psi_lasso_file <- file.path(prefix, paste0(prefix, "_graph_B_lassores_full.csv"))

Omega_file <- file.path(prefix, paste0(prefix, "_graph_Omega_full.csv"))
Omega_lasso_file <- file.path(prefix, paste0(prefix, "_graph_Omega_lassores_full.csv"))

if(!file.exists(Psi_file)) stop("B file does not exist")
if(!file.exists(Omega_file)) stop("Omega file does not exist")

Psi_results <- rbind(read.csv(Psi_file), read.csv(Psi_lasso_file))
Omega_results <- rbind(read.csv(Omega_file), read.csv(Omega_lasso_file))


dims <- paste0("p",p,"q",q,"n",n,".png")

png(paste0("Psi_Omega_sen_prec_frob_", dims), width = 11, height = 6, units = "in", res = 400)
par(mar = c(1,1,1.5,1), mgp = c(0, 0.3, 0), mfrow = c(5,6),
    cex.main = 1.5)

model_names <- c("AR(1)","AR(2)","Block","Star","Dense")

for(i in 1:5){
tmp <- get_result_list(Psi_results[Psi_results$mod==(i+(i==5)),], 
                                 Omega_results[Omega_results$mod==(i+(i==5)),])

#png("Psi_Omega_sen_prec.png", width = 6, height = 6, units = "in", res = 400)
#par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,2))
#
cat("\n",model_names[i],"\n")
get_result_table(tmp)
plot_result_list(tmp, model_names[i], i==1)

}

dev.off()
