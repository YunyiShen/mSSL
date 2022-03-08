library(mSSL)
gut <- CARlasso::mgp154
soil <- CARlasso::mgp2592

gut_Y <- as.matrix( gut[,2:16])
gut_Y <- gut_Y + (gut_Y==0)
gut_Y <- apply( gut_Y[,1:14], 2, `/`, gut_Y[,15] )
gut_Y_transform <- log(gut_Y)
gut_X <- as.matrix( model.matrix(~., gut[,17:21]))[,-1]
gut_X <- apply(gut_X,2,function(w){(w-mean(w))/sd(w)})

gut_cgdperes <- mSSL(gut_Y_transform, gut_X, verbose = T, eps = 1e-4, diag_penalty = FALSE, max_iter = 10000)
gut_mdperes <- mSSL(gut_Y_transform, gut_X, cg = FALSE, verbose = T)

library(caret)
set.seed(42)
folds_gut <- createFolds(1:nrow(gut_Y_transform), k = 10)

tenfold_gut <- lapply(folds_gut, function(idx, X,Y){
  x_train <- X[-idx,]
  y_train <- Y[-idx,]
  x_test <- X[idx,]
  y_test <- Y[idx,]
  gut_cgdperes <- mSSL(y_train, x_train, verbose = F, eps = 1e-4, diag_penalty = FALSE, max_iter = 10000)
  gut_mdperes <- mSSL(y_train, x_train, cg = FALSE, verbose = F)
  
  res <- data.frame(matrix(NA,1,2))
  colnames(res) <- c("mSSL","cgSSL")
  
  res[1,1] <- mean((y_test - x_test %*% gut_mdperes$B - rep(1,nrow(x_test)) %*% t(gut_mdperes$alpha) )^2)
  Sigma <- solve(gut_cgdperes$Omega )
  res[1,2] <- mean((y_test - x_test %*% gut_cgdperes$B %*% Sigma  - rep(1,nrow(x_test)) %*% t(gut_cgdperes$alpha) %*% Sigma )^2)
  return(res)
},gut_X, gut_Y_transform )




soil_Y <- as.matrix( soil[,2:19])
soil_Y <- soil_Y + (soil_Y==0)
soil_Y <- apply( soil_Y[,1:17], 2, `/`, soil_Y[,18] )
soil_Y_transform <- log(soil_Y)
soil_X <- as.matrix(model.matrix(~., soil[,20:26]))[,-c(1,2,3,5,7)]
soil_X <- apply(soil_X, 2, function(w){(w-mean(w))/sd(w)})

soil_cgdperes <- mSSL(soil_Y_transform, soil_X, verbose = T, s_max_condition = 5000, diag_penalty = F)
soil_mdperes <- mSSL(soil_Y_transform, soil_X, cg = FALSE, verbose = T, s_max_condition = 5000, diag_penalty = F)


set.seed(42)
folds_soil <- createFolds(1:nrow(soil_Y_transform), k = 10)

tenfold_soil <- lapply(folds_soil, function(idx, X,Y){
  x_train <- X[-idx,]
  y_train <- Y[-idx,]
  x_test <- X[idx,]
  y_test <- Y[idx,]
  soil_cgdperes <- mSSL(y_train, x_train, verbose = F, eps = 1e-4, diag_penalty = FALSE, max_iter = 10000)
  soil_mdperes <- mSSL(y_train, x_train, cg = FALSE, verbose = F)
  
  res <- data.frame(matrix(NA,1,2))
  colnames(res) <- c("mSSL","cgSSL")
  
  res[1,1] <- mean((y_test - x_test %*% soil_mdperes$B - rep(1,nrow(x_test)) %*% t(soil_mdperes$alpha) )^2)
  Sigma <- solve(soil_cgdperes$Omega )
  res[1,2] <- mean((y_test - x_test %*% soil_cgdperes$B %*% Sigma  - rep(1,nrow(x_test)) %*% t(soil_cgdperes$alpha) %*% Sigma )^2)
  return(res)
},soil_X, soil_Y_transform )

tenfold_gut_flat <- Reduce(rbind, tenfold_gut)
tenfold_soil_flat <- Reduce(rbind, tenfold_soil)

save.image("./real_data/Micorbiome/res_cv.RData")
