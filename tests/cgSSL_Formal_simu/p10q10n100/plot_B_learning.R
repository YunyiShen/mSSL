library(ggplot2)

modname <- c("AR1","AR2","Block","Star","Full")

all_B_q10 <- read.csv("./p10q10n100_graph_B_full.csv")
all_B_q10$mod[all_B_q10$mod==6] <- 5

all_B_q10$beta.sparsity <- factor( 1-all_B_q10$s,levels = c(0.8))
all_B_q10$p <- paste0(all_B_q10$p, " predictors")
all_B_q10 <- within(all_B_q10, p<-factor(p, levels=c( "10 predictors")))
all_B_q10$mod <- modname[all_B_q10$mod]
all_B_q10$mod <- factor(all_B_q10$mod, levels = modname)

B_learning_q10 <- all_B_q10
B_learning_q10[is.na(B_learning_q10)] <- 0
#B_learning_q10$beta.sparsity <- B_learning_q10$s

## New plot
B_learning_q10 <- within(B_learning_q10, algo<-factor(algo, levels=c("cgSSL-dpe", "cgSSL-dcpe", "mSSL-dpe","mSSL-dcpe", "CAR-A", "CAR" )))
B_MCC <- ggplot(data = B_learning_q10,aes(x=algo,y = MCC)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot() + 
  facet_grid(~mod) + 
  ylab("MCC on B") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

B_MCC
ggsave("./B_MCC.jpg",B_MCC,width = 8,height = 5,unit = "in")


B_Sensitivity <- ggplot(data = B_learning_q10[B_learning_q10$mod!="model 5",],aes(x=algo,y = SEN)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot() + 
  facet_grid(~mod) + 
  ylab("Sensitivity on B") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

B_Sensitivity


B_Specificity <- ggplot(data = B_learning_q10[B_learning_q10$mod!="model 5",],aes(x=algo,y = SPE)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot() + 
  facet_grid(~mod) + 
  ylab("Specificity on B") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

B_Specificity


B_FROB <- ggplot(data = B_learning_q10,aes(x=algo,y = log(FROB))) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot() + 
  facet_grid(~mod) + 
  ylab("log FROB on B") + 
  xlab("") + theme_bw() +
  theme(legend.position="top") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

B_FROB
