library(dplyr)
library(tidyr)
library(ggplot2)

ESP1 <- function(a,n){
  return(1+a/n)
}
ESP2 <- function(a,n){
  return(1+a/n^(3/4))
}
ESP1_M <- function(a,n){
  X <- c()
  for (i in 1:length(n)) {
    X[i] <- rnorm(1,0,00.5/n[i])
  }
  return(1+a/n+X)
}
ESP2_M <- function(a,n){
  X <- c()
  for (i in 1:length(n)) {
    X[i] <- rnorm(1,0,00.5/n[i])
  }
  return(1+a/n^(3/4)+X)
}

N.grid <-seq(10,100,1)

d <- data.frame("Samplesize"=N.grid,"ESPE1_2"=ESP1(2,N.grid),"ESPE1_3"=ESP1(3,N.grid),"ASPE1_2"=ESP1_M(2,N.grid),"ASPE1_3"=ESP1_M(3,N.grid))

ggplot(d,aes(x=Samplesize)) +
  geom_line(aes(y=ESPE1_2,color = "M*"),size =1,linetype=5) +
  geom_line(aes(y=ASPE1_2,color = "M*"),size =1.2) +
  geom_line(aes(y=ESPE1_3,color = "M*+1"),size =1,linetype=5) +
  geom_line(aes(y=ASPE1_3,color = "M*+1"),size =1.2) +
  ylab("Expected Squared Prediction Error") + ggtitle("Asymptotic Properties CV(1)") +
  theme_bw() +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20,hjust = 0.5,face="bold")) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(plot.title = element_text(margin = margin(t = 0, r = 0 , b = 20, l = 0))) +
  theme(legend.title=element_blank(),legend.text = element_text(size = 13)) +
  theme(legend.position = c(0.85,0.8))

d <-data.frame("Samplesize"=N.grid,"ESPE2_2"=ESP2(2,N.grid),"ESPE2_3"=ESP2(3,N.grid),"ASPE2_2"=ESP2_M(2,N.grid),"ASPE2_3"=ESP2_M(3,N.grid))

ggplot(d,aes(x=Samplesize)) +
  geom_line(aes(y=ESPE2_2,color = "M*"),size =1,linetype=5) +
  geom_line(aes(y=ASPE2_2,color = "M*"),size =1.2) +
  geom_line(aes(y=ESPE2_3,color = "M*+1"),size =1,linetype=5) +
  geom_line(aes(y=ASPE2_3,color = "M*+1"),size =1.2) +
  ylab("Expected Squared Prediction Error") + ggtitle(expression(bold(paste("Asymptotic Properties ",CV(n[v]))))) +
  theme_bw() +
  theme(text = element_text(size=15),
        plot.title = element_text(size = 20,hjust = 0.5,face="bold")) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  theme(plot.title = element_text(margin = margin(t = 0, r = 0 , b = 20, l = 0))) +
  theme(legend.title=element_blank(),legend.text = element_text(size = 13)) +
  theme(legend.position = c(0.85,0.8))






plot(N.grid,ESP1(2,N.grid),type = "l")
lines(N.grid,ESP1_M(2,N.grid),lty =2)
lines(N.grid,ESP1(3,N.grid),col="red")
lines(N.grid,ESP1_M(3,N.grid),lty = 2,col = "red")

plot(N.grid,ESP1(2,N.grid),type = "l",lty=2)
lines(N.grid,ESP1_M(2,N.grid))
lines(N.grid,ESP1(3,N.grid),col="red",lty =2)
lines(N.grid,ESP1_M(3,N.grid),col = "red")


plot(N.grid,ESP2(2,N.grid),type = "l",lty=2)
lines(N.grid,ESP2_M(2,N.grid))
lines(N.grid,ESP2(3,N.grid),col="red",lty =2)
lines(N.grid,ESP2_M(3,N.grid),col = "red")


