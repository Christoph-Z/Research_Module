library(ggplot2)

n <- 100
p <- 5
p.True <- 2
beta <- c(4,7,0,0,0)
DataGen_X <- function(n){
  x_1 <- rep(1,n)         
  x_2 <- rnorm(n,0,1)
  x_3 <- rnorm(n,0,1)
  x_4 <- rnorm(n,0,1)
  x_5 <- rnorm(n,4,1)
  return(X<- cbind(x_1,x_2,x_3,x_4,x_5))
}
X <- DataGen_X(n)
sigma <- 1
y <- X%*%beta+rnorm(n,0,sigma)

##--------------------------------------------------------------------
##ESPE for CV(1) and CV(n_v) for Models in Category II

ESPE_CV1 <- function(d,n){
  return((1+d/n)*sigma^2)
}
ESPE_CVn.v <- function(d,n){
  return(1+d/n^(3/4))
}

##--------------------------------------------------------------------
##ESPE for CV(1) and CV(n_v) for Models in Category II

CV <- function(n,y,X,Method = "CV1"){
  if(Method =="CV1"){
    n_v <- 1
    b <- n
    train <- combn(seq(1,n,1),n-n_v)
  }else{
    b <- Method
    n_v <- n-floor(n^{3/4})
    train <- matrix(ncol = b, nrow = n-n_v)
    for (i in 1:b) {
      train[,i] <- sample(seq(1,n,1),n-n_v,replace = FALSE)
    }
  }
  
  Pred.Error <- c()
  for (i in 1:length(train[1,])) {
    #To make the Code a bit shorter
    train.i <- train[,i]                 
    X.train <- X[train.i,]    
    #Prediction Error for a given alpha and a given subset
    Error <- y[-train.i]-X[-train.i,]%*%solve(t(X.train)%*%X.train)%*%t(X.train)%*%y[train.i]
    Pred.Error[i] <- t(Error)%*%Error
  }
  MeanPred.Error <- 1/(b*n_v)*sum(Pred.Error)
  return(MeanPred.Error)
}


N.grid <-seq(20,100,1)

CV1_d2 <- c()
for (i in 1:length(N.grid)) {
  CV1_d2[i] <- CV(N.grid[i],y[1:N.grid[i]],X[1:N.grid[i],1:2])
}

CV1_d3 <- c()
for (i in 1:length(N.grid)) {
  CV1_d3[i] <- CV(N.grid[i],y[1:N.grid[i]],X[1:N.grid[i],1:3])
}

MCCV1_d2 <- c()
for (i in 1:length(N.grid)) {
  MCCV1_d2[i] <- CV(N.grid[i],y[1:N.grid[i]],X[1:N.grid[i],1:2],Method = N.grid[i]*5)
}

MCCV1_d3 <- c()
for (i in 1:length(N.grid)) {
  MCCV1_d3[i] <- CV(N.grid[i],y[1:N.grid[i]],X[1:N.grid[i],1:3],Method = N.grid[i]*5)
}

##--------------------------------------------------------------------
##Plot for CV1

d <- data.frame("Samplesize"=N.grid,"ESPE1_2"=ESPE_CV1(2,N.grid),"ESPE1_3"=ESPE_CV1(3,N.grid),"ASPE1_2"=CV1_d2,"ASPE1_3"=CV1_d3)

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

##--------------------------------------------------------------------
##Plot CV(n_v)


d <- data.frame("Samplesize"=N.grid,"ESPE1_2"=ESPE_CVn.v(2,N.grid),"ESPE1_3"=ESPE_CVn.v(3,N.grid),"ASPE1_2"=MCCV1_d2,"ASPE1_3"=MCCV1_d3)

ggplot(d,aes(x=Samplesize)) +
  geom_line(aes(y=ESPE1_2,color = "M*"),size =1,linetype=5) +
  geom_line(aes(y=ASPE1_2,color = "M*"),size =1.2) +
  geom_line(aes(y=ESPE1_3,color = "M*+1"),size =1,linetype=5) +
  geom_line(aes(y=ASPE1_3,color = "M*+1"),size =1.2) +
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