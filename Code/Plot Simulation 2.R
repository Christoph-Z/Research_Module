Crit <- rep(c("CV","AIC","MCCV","BIC"),4)
Model <- c(rep("M*",4),rep("M*+1",4),rep("M*+2",4),rep("M*+3",4))
Prob <- t(matrix(c(0.5980,0.336,0.0590,0.0070,0.6015,0.3320,0.0590,0.0075,0.9510,0.048988,0,0,0.9690,0.305,0.0005,0),4,4))
d <-data.frame("C"=Crit,"M"=Model,"P"=c(0.5980,0.6015,0.9510,0.9690,0.3366,0.3320,0.048988,0.0305,0.0590,0.0590,0,0.0005,0.007,0.007,0,0),"COL"=rep("blue",16),"Prozent"=c("60%","60%","95%","97%","34%","33%","4%","3%","6%","6%","0%","0%","0.7%","0.8%","0%","0%"))

ggplot(d,aes(x=C,y=M,size=P,label=P)) +
  geom_point(alpha=0.8,stat="identity",color ='skyblue3') +
  xlab("Selection Method") + ylab("Model Dimensionality") + ggtitle("Probabilties of Predicting Cat. II Model") +
  scale_size_continuous(range = c(1, 35)) +
  theme_bw() +
  guides(colour = FALSE,size = FALSE) +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20,hjust = 0.5,face="bold")) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(plot.title = element_text(margin = margin(t = 0, r = 0 , b = 20, l = 0))) +
  geom_text(size =5,aes(label=ifelse(P>0.04,as.character(Prozent),''))) +
  geom_text(size =5,aes(label=ifelse(P<0.04,as.character(Prozent),'')),hjust=-0.5,vjust=0.5)




