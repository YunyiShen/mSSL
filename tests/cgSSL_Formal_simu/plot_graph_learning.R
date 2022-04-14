library(ggplot2)
modname <- c("AR1","AR2","Block","Star","Full")
all_Omega_p10 <- read.csv("./p100q30n400/p100q30n400_graph_Omega_full.csv")
all_Omega_p10$mod[all_Omega_p10$mod==6] <- 5

all_Omega_p10$beta.sparsity <- factor( 1-all_Omega_p10$s,levels = c(0.8))
all_Omega_p10$p <- paste0(all_Omega_p10$p, " predictors")
all_Omega_p10 <- within(all_Omega_p10, p<-factor(p, levels=c( "100 predictors")))
all_Omega_p10$mod <- modname[all_Omega_p10$mod]
all_Omega_p10$mod <- factor(all_Omega_p10$mod, levels = modname)

Omega_learning_p10 <- all_Omega_p10
Omega_learning_p10[is.na(Omega_learning_p10)] <- 0
#Omega_learning_p10$beta.sparsity <- Omega_learning_p10$s

## New plot
Omega_learning_p10 <- within(Omega_learning_p10, algo<-factor(algo, levels=rev(c("cgSSL-dpe", "cgSSL-dcpe", "mSSL-dpe","mSSL-dcpe", "CAR-A", "CAR" ))))
Graph_MCC <-  ggplot(data = Omega_learning_p10[Omega_learning_p10$mod!="Full",],aes(x=algo,y = F1)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~mod) + 
  ylab("F1 on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        #axis.text.y = element_text(angle=180),
        plot.margin = margin(.15, .15, .15, .15, "cm")) + 
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + coord_flip()

Graph_MCC
ggsave("./p100q30n400/Omega_F1.pdf",Graph_MCC,width = 9,height = 5,unit = "in")


Graph_Sensitivity <- ggplot(data = Omega_learning_p10[Omega_learning_p10$mod!="Full",],aes(x=algo,y = SEN)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~mod) + 
  ylab("Sensitivity on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        #axis.text.y = element_text(angle=180),
        plot.margin = margin(.15, .15, .15, .15, "cm")) + 
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + coord_flip()

Graph_Sensitivity


Graph_Specificity <- ggplot(data = Omega_learning_p10[Omega_learning_p10$mod!="Full",],aes(x=algo,y = SPE)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~mod) + 
  ylab("Specificity on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        #axis.text.y = element_text(angle=180),
        plot.margin = margin(.15, .15, .15, .15, "cm")) + 
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + coord_flip()

Graph_Specificity


Graph_FROB <- ggplot(data = Omega_learning_p10[Omega_learning_p10$mod!="Full",],aes(x=algo,y = log(FROB))) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~mod) + 
  ylab("log FROB loss on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        #axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text.x = element_text(angle=45,hjust = 1,vjust=1),
        #axis.text.y = element_text(angle=180),
        plot.margin = margin(.15, .15, .15, .15, "cm")) + 
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + coord_flip()

Graph_FROB
ggsave("./p100q30n400/Omega_FROB.pdf",Graph_FROB,width = 9,height = 5,unit = "in")

