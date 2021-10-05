gut <- CARlasso::mgp154
soil <- CARlasso::mgp2592

gut_Y <- as.matrix( gut[,2:16])
gut_Y <- apply( gut_Y[,1:14], 2, `/`, gut_Y[,15] )
gut_Y_transform <- asinh(gut_Y)
gut_X <- as.matrix( model.matrix(~.-1, gut[,17:21]))

gut_cgdperes <- mSSL(gut_Y_transform, gut_X, verbose = T)
gut_mdperes <- mSSL(gut_Y_transform, gut_X, cg = FALSE, verbose = T)

soil_Y <- as.matrix( soil[,2:19])
soil_Y <- apply( soil_Y[,1:17], 2, `/`, soil_Y[,18] )
soil_Y_transform <- asinh(soil_Y)
soil_X <- as.matrix(model.matrix(~.-1, soil[,19:26]))

soil_cgdperes <- mSSL(soil_Y_transform, soil_X, verbose = T, s_max_condition = 5000)
soil_mdperes <- mSSL(soil_Y_transform, soil_X, cg = FALSE, verbose = T)


