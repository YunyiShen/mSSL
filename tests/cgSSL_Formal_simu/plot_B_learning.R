library(ggplot2)

modname <- c("AR1","AR2","Block","Star","Full")

all_B_q10 <- rbind(read.csv("./p10q10n100/p10q10n100_graph_B_full.csv"), 
                   read.csv("./p10q10n100/p10q10n100_graph_B_lassores_full.csv"))
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
B_learning_q10 <- within(B_learning_q10, algo<-factor(algo, levels= rev( c("cgSSL-dpe", "cgSSL-dcpe", "mSSL-dpe","mSSL-dcpe", "CAR-A", "CAR", "LASSO" ))))
B_MCC <- ggplot(data = B_learning_q10,aes(x=algo,y = F1)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~mod) + 
  ylab("F1 on B") + 
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

B_MCC 
ggsave("./p10q10n100/B_F1.pdf",B_MCC,width = 10,height = 5,unit = "in")


B_Sensitivity <-  ggplot(data = B_learning_q10,aes(x=algo,y = SEN)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~mod) + 
  ylab("Sensitivity on B") + 
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

B_Sensitivity


B_Specificity <-  ggplot(data = B_learning_q10,aes(x=algo,y = SPE)) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~mod) + 
  ylab("Specificity on B") + 
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

B_Specificity


B_FROB <-  ggplot(data = B_learning_q10,aes(x=algo,y = log(FROB))) + 
  geom_point( alpha=0.1, size=1)+
  geom_boxplot(linetype = "dashed", outlier.shape = 1) + 
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) +
  stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width = 0.5)+
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width = 0.5)+
  facet_grid(~mod) + 
  ylab("log FROB loss on B") + 
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

B_FROB
ggsave("./p100q30n400/B_FROB.pdf",B_FROB,width = 10,height = 5,unit = "in")

