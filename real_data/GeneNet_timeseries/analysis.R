library(mSSL)
nine_genes <- read.csv("./real_data/GeneNet_timeseries/nine_gene.csv", row.names = 1)
nine_genes <- apply(nine_genes, 2, function(w){w-mean(w)})[1:6,]
dpe_res <- mSSL::cgVARSSL(nine_genes,verbose = T, s_max_condition = 5000)
save.image("./real_data/GeneNet_timeseries/cgVARres.RData")
