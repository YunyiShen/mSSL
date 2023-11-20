library(igraph)
library(GGally)
library(Matrix)
library(ggplot2)
library(ggraph)
source("./real_data/misc.R")
# Some resources on network plots: http://mr.schochastics.net/netVizR.html


plot_graph <- function(graph){
  dd <- c(rep(0,6),0.1,rep(0,12),-0.1,rep(0,5))
  p <- ggraph(graph,layout = "circle")+
    geom_edge_link(
      aes(color = edge.color), width = 1.5,alpha = 1, 
      arrow = arrow(angle = 15, length = unit(0.15, "inches"), type = "closed")) + 
    scale_edge_color_manual(values = (c("black","red")))+
    geom_node_point(mapping = aes(shape = shape,
                                  stroke = 1.5),
                    size = 5,
                    col = "black",
                    fill= "lightgray", alpha = 1) + 
    scale_shape_manual(values = c(21,24)) + 
    coord_fixed(clip = 'off')+
    guides(width = guide_legend(order = 1),
           #color = guide_legend(order=3),
           size = guide_legend(order=2),
           shape = "none", #guide_legend(order=1),
           edge_color = "none" #guide_legend(order = 3)
           #alpha = FALSE, linetype = FALSE
    ) 
  p <- p+
    geom_node_text(aes(label = name),nudge_x = p$data$x * .38, 
                   nudge_y = p$data$y * .2+dd, family = "")+  #repel = T,check_overlap = T)+
    theme_graph(base_family = 'Helvetica')
}



# for human
load("./real_data/Micorbiome/res.RData")
load("./real_data/Micorbiome/res_bb.RData")

signif_psi <- matrix(FALSE, 11,14)
signif_omega <- matrix(FALSE, 14,14)
for(j in 1:14){
  for(i in 1:11){
    tmp <- gut_cg_BB$samples$Psi[i,j,]
    signif_psi[i,j] <- (quantile(tmp, .025, na.rm = T) * 
                          quantile(tmp, .975, na.rm = T)) > 0
  }
  for(i in 1:14){
    tmp <- gut_cg_BB$samples$Omega[i,j,]
    signif_omega[i,j] <- (quantile(tmp, .025, na.rm = T) * 
                            quantile(tmp, .975, na.rm = T)) > 0
  }
}


col_pn <- c("lightgray","red")

## vertices list
vertices_df <- data.frame(id=c(paste0("M",1:(ncol(gut_Y))),
                               paste0("E",1:ncol(gut_X))),
                          group = c(rep("microbe",(ncol(gut_Y))),
                                    rep("env",ncol(gut_X))),
                          shape = c(rep("circle",(ncol(gut_Y))),
                                    rep("square",ncol(gut_X))))
## edges between microbes
ind_mat_micro <- expand.grid(from = 1:(ncol(gut_Y)),to = 1:(ncol(gut_Y)))
#ind_mat_micro <- ind_mat_micro[ind_mat_micro$from>ind_mat_micro$to,]
ind_mat_micro$weight <- sapply(1:nrow(ind_mat_micro),function(i,indmat,mat){
  mat[indmat$from[i],indmat$to[i]]
},ind_mat_micro,gut_cgdperes$Omega) # get the weight

### significance via BB
ind_mat_micro$signif <- sapply(1:nrow(ind_mat_micro),function(i,indmat,mat){
  mat[indmat$from[i],indmat$to[i]]
},ind_mat_micro,signif_omega) # get signif

## edges between env and microbes
ind_mat_env <- expand.grid(from = 1:ncol(gut_X),to = 1:(ncol(gut_Y)))
ind_mat_env$weight <- sapply(1:nrow(ind_mat_env),function(i,indmat,mat){
  mat[indmat$from[i],indmat$to[i]]
},ind_mat_env,gut_cgdperes$B)

ind_mat_env$signif <- sapply(1:nrow(ind_mat_env),function(i,indmat,mat){
  mat[indmat$from[i],indmat$to[i]]
},ind_mat_env,signif_psi) # get signif


# get directions and type of edges (by color)
ind_mat_env$from <- paste0("E",ind_mat_env$from)
ind_mat_env$to <- paste0("M",ind_mat_env$to)
ind_mat_env$col <- col_pn[2]
ind_mat_micro$from <- paste0("M",ind_mat_micro$from)
ind_mat_micro$to <- paste0("M",ind_mat_micro$to)
ind_mat_micro$col <- col_pn[1]

# combined edges
edge_df <- rbind(ind_mat_micro,ind_mat_env)
edge_df <- edge_df[edge_df$weight!=0 | edge_df$signif,]

full_graph <- graph.data.frame(edge_df,vertices_df,directed=T)


#col_ER <- c("orange","darkgreen")
shape_ER <- c("square","circle")
type <- c("predictors","microbe")
direction <- c("negative","positive")

## some edge settings
E(full_graph)$edge.color <- E(full_graph)$col
E(full_graph)$non0 <- abs( E(full_graph)$weight)>0
E(full_graph)$linetypes <- c("g","a")[E(full_graph)$signif + 1]
V(full_graph)$name <- c(colnames(gut_Y),colnames(gut_X))
V(full_graph)$name[16] <- "Male"
V(full_graph)$name[17] <- "DayHospital"
V(full_graph)$name[18] <- "Long-term"
V(full_graph)$name[19] <- "Rehab"
V(full_graph)$name[20] <- "Diet1"
V(full_graph)$name[21] <- "Diet2"
V(full_graph)$name[22] <- "Diet3"
V(full_graph)$name[23] <- "Diet4"
V(full_graph)$name[24] <- "DietPEG"


p <- plot_graph(full_graph)
p
ggsave("./real_data/Micorbiome/humangut2_BA.pdf",p,width = 7,height = 7, scale = .75)

siggraph <- full_graph
siggraph <- delete_edges(siggraph, E(siggraph)[!E(siggraph)$signif] )
p2 <- plot_graph(siggraph)
ggsave("./real_data/Micorbiome/humangut_signig2_BA.pdf",p2,width = 7,height = 7, scale = .75)
