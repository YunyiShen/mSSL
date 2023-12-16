library(igraph)
library(GGally)
source("./real_data/misc.R")
# Some resources on network plots: http://mr.schochastics.net/netVizR.html

# for human
load("./real_data/GeneNet_timeseries/cgVARres.RData")


among_spp <- graph_from_adjacency_matrix(dpe_res$Omega,mode = "directed",weighted = T,diag = F)
linear_reg_graph <- graph_from_adjacency_matrix(t(dpe_res$B),weighted = T, mode = "directed")

col_pn <- c("lightblue","pink")
l <-layout.circle(among_spp)#, repulserad=vcount(among_spp)^3,area=vcount(among_spp)^2.4)

pdf("./real_data/GeneNet_timeseries/res.pdf", width = 4, height = 4)
par(mar=c(0,0,0,0)+.1)
plot(among_spp,edge.arrow.size=.1,
     vertex.label=colnames(nine_genes),
     #vertex.size = (3* alpha_centrality(among_spp)),
     #vertex.label=1:35,
     layout = layout.circle,
     vertex.color = c("red","red","yellow","yellow","red","yellow","yellow","red","red"),
     edge.curved=0)
plot(linear_reg_graph,edge.arrow.size=.5,
     vertex.label=colnames(nine_genes),
     add = T,
     #vertex.size = (3* alpha_centrality(among_spp)),
     #vertex.label=1:35,
     layout = layout.circle,
     vertex.color = c("red","red","yellow","yellow","red","yellow","yellow","red","red"),
     #color = c("red","red","yellow","yellow","red","yellow","yellow","red","red"),
     edge.color = "blue",
     edge.curved=.3)
dev.off()
