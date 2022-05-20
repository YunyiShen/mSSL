p <- 10
q <- 10
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

## WHY DOES Psi_results and Omega_results have different number of rows??? ###

algo <- c("LASSO", "CAR", "CAR-A", "cgSSL-dcpe", "cgSSL-dpe")


  
outcomes <- c("SEN", "PREC")
for(out in outcomes){
  tmp <- data.frame("LASSO" = rep(NA, times = 501),
                    "CAR" = rep(NA, times = 501),
                    "CAR-A" = rep(NA, times = 501),
                    "cgSSL-dcpe" = rep(NA, times = 501),
                    "cgSSL-dpe" = rep(NA, times = 501))
  for(i in 1:5){
    index <- which(Psi_results[,"algo"] == algo[i] & Psi_results[,"mod"]==1)
    #if(algo[i] != "LASSO") tmp[,i] <- Psi_results[index,out]
    #else tmp[1:500,i] <- Psi_results[index, out]
    tmp[1:length(index),i] <- Psi_results[index, out]
  }
  assign(paste0(out, "_psi_results"), tmp)
  
  tmp <- data.frame("LASSO" = rep(NA, times = 501),
                    "CAR" = rep(NA, times = 501),
                    "CAR-A" = rep(NA, times = 501),
                    "cgSSL-dcpe" = rep(NA, times = 501),
                    "cgSSL-dpe" = rep(NA, times = 501))
  for(i in 1:5){
    index <- which(Omega_results[,"algo"] == algo[i] & Omega_results[,"mod"] == 1)
    #if(algo[i] != "LASSO") tmp[,i] <- Omega_results[index,out]
    #else tmp[1:500,i] <- Omega_results[index, out]
    tmp[1:length(index),i] <- Omega_results[index,out]
  }
  assign(paste0(out, "_omega_results"), tmp)
  
}

png("Psi_sen_prec.png", width = 6, height = 3, units = "in", res = 400)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))
boxplot(SEN_psi_results, horizontal = TRUE, main = expression(Psi~"sensitivity"),
        medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
        ylim = c(0,1), yaxt = "n")
text(x = 0.2, y = 1, labels = "cgLASSO")
text(x = 0.2, y = 2, labels = "CAR")
text(x = 0.2, y = 3, labels = "CAR-A")
text(x = 0.2, y = 4, labels = "cgSSL-DCPE")
text(x = 0.2, y = 5, labels = "cgSSL-DPE")
boxplot(PREC_psi_results, horizontal = TRUE, main = expression(Psi~"precision"),
        medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
        ylim = c(0,1), yaxt = "n")
#text(x = 0.2, y = 1, labels = "LASSO")
#text(x = 0.2, y = 2, labels = "CAR")
#text(x = 0.2, y = 3, labels = "CAR-A")
#text(x = 0.2, y = 4, labels = "cgSSL-DCPE")
#text(x = 0.2, y = 5, labels = "cgSSL-DPE")
dev.off()


#png("Psi_Omega_sen_prec.png", width = 6, height = 6, units = "in", res = 400)
#par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,2))
png("Psi_Omega_sen_prec.png", width = 8, height = 2.2, units = "in", res = 400)
par(mar = c(3,1,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,4),
    cex.main = 1.5)

boxplot(SEN_psi_results, horizontal = TRUE, main = expression(Psi~"sensitivity"),
        medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
        ylim = c(0,1), yaxt = "n")
text(x = 0.2, y = 1, labels = "cgLASSO", cex = 1.2)
text(x = 0.2, y = 2, labels = "CAR", cex = 1.2)
text(x = 0.2, y = 3, labels = "CAR-A", cex = 1.2)
text(x = 0.22, y = 4, labels = "cgSSL-DCPE", cex = 1.2)
text(x = 0.2, y = 5, labels = "cgSSL-DPE", cex = 1.2)
boxplot(PREC_psi_results, horizontal = TRUE, main = expression(Psi~"precision"),
        medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
        ylim = c(0,1), yaxt = "n")

boxplot(SEN_omega_results, horizontal = TRUE, main = expression(Omega~"sensitivity"),
        medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
        ylim = c(0,1), yaxt = "n")

boxplot(PREC_omega_results, horizontal = TRUE, main = expression(Omega~"precision"),
        medlwd = 0.5, pch = 16, cex = 0.75, names = NA,
        ylim = c(0,1), yaxt = "n")

dev.off()
