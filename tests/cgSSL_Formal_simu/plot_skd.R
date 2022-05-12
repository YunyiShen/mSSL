p <- 10
q <- 10
n <- 100

prefix <- paste0("p", p, "q", q, "n", n)
B_file <- file.path(prefix, paste0(prefix, "_graph_B_full.csv"))
B_lasso_file <- file.path(prefix, paste0(prefix, "_graph_B_lassores_full.csv"))

Omega_file <- file.path(prefix, paste0(prefix, "_graph_Omega_full.csv"))

if(!file.exists(B_file)) stop("B file does not exist")
if(!file.exists(Omega_file)) stop("Omega file does not exist")

B_results <- read.csv(B_file)
B_lasso_results <- read.csv(B_lasso_file)


B_results <- rbind(B_results, B_lasso_results)
for(mod in unique(B_results[,"mod"])){
  index <- which(B_results[,"mod"] == mod)
  assign(paste0("B_results_", mod), 
         B_results[index,])
}

algo <- c("CAR", "CAR-A", "cgSSL-dcpe", "cgSSL-dpe", "LASSO")
tmp <- data.frame("CAR" = rep(NA, times = 500),
                  "CAR-A" = rep(NA, times = 500),
                  "cgSSL-dcpe" = rep(NA, times = 500),
                  "cgSSL-dpe" = rep(NA, times = 500),
                  "LASSO" = rep(NA, times = 500))
for(a in algo){
  print(a)
  index <- which(B_results[,"algo"] == a)
  tmp_results <- B_results[index,]
  tmp[,a] <- tmp_results[1:500, "F1"]
}

boxplot(tmp[,algo])
