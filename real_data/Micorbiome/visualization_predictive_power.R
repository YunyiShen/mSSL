library(ggplot2)
library(reshape2)
ggplot(data = reshape2::melt( tenfold_gut_flat), aes(y=value,x=variable)) + 
  geom_boxplot() + 
  geom_point() + 
  xlab("") + ylab("relative MSE of 10-fold CV") + 
  theme_classic()

ggsave("predictive_power_gut.pdf", width = 6, height = 3, scale = 0.8)
