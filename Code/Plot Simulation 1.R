library(dplyr)
library(tidyr)

MCV <- MCV[1:50]
AIC <- AIC[1:50]
CV1 <- CV1[1:50]
BIC <- BIC[1:50]
N <- N[1:50]

plot(N,BIC,type = "l",lwd = 2 ,ylim = c(0.62,1),xlim = c(15,113), xlab = "Samplesize", ylab = "Prob. of choosing a Cat II Model", main = "Probability of Cat II Model")
lines(N,AIC,type = "l",lwd = 2,col = "red")
lines(N,CV1,type = "l",lwd = 2, col = "green")
lines(N,MCV,type = "l",lwd = 2,col = "blue")
axis(1,at = seq(10,130,20))
axis(2,at = c(0.5,0.6,0.7,0.8,0.9,1,1.5))
par(lwd=1)
legend(13,1,c("BIC","AIC","CV(1)","MCCV"),col = c("black","red","green","blue"),lty = 1,lwd = 2,cex = 1)


d <- data.frame("Samplesize"=N,"MCCV"=MCV,"AIC"=AIC,"CV1"=CV1,"BIC"=BIC)
d1 <-d %>% gather(key,value,MCCV,AIC,BIC,CV1)

ggplot(d1,aes(x=Samplesize,y=value,colour=key)) +
  geom_line(size=1.2) +
  coord_cartesian(ylim = c(0.6,1), xlim = c(13, 113)) +
  ylab("Prob of choosing a Cat. II Model") + ggtitle("Probabiltie of Cat. II Model") +
  theme_bw() +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20,hjust = 0.5,face="bold")) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(plot.title = element_text(margin = margin(t = 0, r = 0 , b = 20, l = 0))) +
  theme(legend.title=element_blank(),legend.text = element_text(size = 13)) +
  theme(legend.position = c(0.1,0.8))
