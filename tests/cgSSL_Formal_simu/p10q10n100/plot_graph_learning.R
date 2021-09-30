library(ggplot2)
modname <- c("AR1","AR2","Block","Star","Full")
all_Omega_p10 <- read.csv("./p10q10n100graph_Omega_full.csv")
all_Omega_p10$mod[all_Omega_p10$mod==6] <- 5

all_Omega_p10$beta.sparsity <- factor( 1-all_Omega_p10$s,levels = c(0.8))
all_Omega_p10$p <- paste0(all_Omega_p10$p, " predictors")
all_Omega_p10 <- within(all_Omega_p10, p<-factor(p, levels=c( "10 predictors")))
all_Omega_p10$mod <- modname[all_Omega_p10$mod]
all_Omega_p10$mod <- factor(all_Omega_p10$mod, levels = modname)

Omega_learning_p10 <- all_Omega_p10
Omega_learning_p10[is.na(Omega_learning_p10)] <- 0
#Omega_learning_p10$beta.sparsity <- Omega_learning_p10$s

## New plot
Omega_learning_p10 <- within(Omega_learning_p10, algo<-factor(algo, levels=c("cgSSL-dpe", "cgSSL-dcpe", "mSSL-dpe","mSSL-dcpe", "CAR-A", "CAR" )))
Graph_MCC <- ggplot(data = Omega_learning_p10[Omega_learning_p10$mod!="Full",],aes(x=algo,y = MCC)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot() + 
  facet_grid(~mod) + 
  ylab("MCC on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Graph_MCC
ggsave("Omega_MCC.jpg",Graph_MCC,width = 8,height = 5,unit = "in")


Graph_Sensitivity <- ggplot(data = Omega_learning_p10[Omega_learning_p10$mod!="Full",],aes(x=algo,y = SEN)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot() + 
  facet_grid(~mod) + 
  ylab("Sensitivity on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Graph_Sensitivity


Graph_Specificity <- ggplot(data = Omega_learning_p10[Omega_learning_p10$mod!="Full",],aes(x=algo,y = SPE)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot() + 
  facet_grid(~mod) + 
  ylab("Specificity on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Graph_Specificity


Graph_FROB <- ggplot(data = Omega_learning_p10[Omega_learning_p10$mod!="Full",],aes(x=algo,y = log(FROB))) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot() + 
  facet_grid(~mod) + 
  ylab("log FROB on Omega") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

Graph_FROB
